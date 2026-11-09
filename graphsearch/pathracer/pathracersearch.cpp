// Copyright 2025 Anton Korobeynikov

// This file is part of Bandage-NG

// Bandage is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// Bandage is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with Bandage.  If not, see <http://www.gnu.org/licenses/>.


#include "pathracersearch.h"

#include "graph/debruijnnode.h"
#include "graphsearch/query.h"
#include "program/globals.h"
#include "program/settings.h"

#include "graph/assemblygraph.h"
#include "graph/gfawriter.h"
#include "io/fileutils.h"

#include <QDir>
#include <QFile>
#include <QMap>
#include <QProcess>
#include <QRegularExpression>
#include <QTemporaryFile>
#include <QTextStream>
#include <QVector>

using namespace search;

PathRacerSearch::PathRacerSearch(const QDir &workDir, QObject *parent)
        : GraphSearch(workDir, parent) {}

bool PathRacerSearch::findTools() {
    if (!findProgram("pathracer", g_settings->pathracerCommand, &m_pathracerCommand,
                     "PathRacer (from SPAdes)"))
        return false;

    return true;
}

// PathRacer takes the assembly graph itself (in GFA format) as its input,
// so "building the database" simply means exporting the graph to a GFA file.
// PathRacer performs its own frame-shift aware alignment of HMMs to the
// graph, so no six-frame translation is needed here.
QString PathRacerSearch::buildDatabase(const AssemblyGraph &graph, bool includePaths) {
    DbBuildFinishedRAII watcher(this);

    m_lastError = "";
    if (!findTools())
        return m_lastError;

    bool atLeastOneSequence = false;
    for (const auto *node : graph.m_deBruijnGraphNodes) {
        if (!node->sequenceIsMissing()) {
            atLeastOneSequence = true;
            break;
        }
    }

    if (!atLeastOneSequence)
        return (m_lastError = "Cannot build the PathRacer database as this graph contains no sequences");

    // Note: includePaths is implicit, as the exported GFA always contains
    // the graph paths and PathRacer uses them as scaffolds.
    if (!gfa::saveEntireGraph(temporaryDir().filePath("graph.gfa"), graph))
        return (m_lastError = "Failed to export the graph in GFA format");

    return m_lastError;
}

QString PathRacerSearch::doSearch(QString extraParameters) {
    return doSearch(queries(), extraParameters);
}

static void writeQueryFile(QFile *file,
                           const Queries &queries) {
    QTextStream out(file);
    for (const auto *query: queries.queries()) {
        out << query->getAuxData()
            << "//\n";
    }
}

namespace {

// A single alignment path as reported by PathRacer in the P-lines of the
// annotated graph.  The P-line name looks like this:
// <query name>_<path index>_score_<score>
struct PathRacerPath {
    QString queryName;
    int idx = -1;
    double bitScore = 0.0;
    QString pathString; // e.g. "1+,2+" or "5-,4-"
};

// A single alignment as reported by PathRacer in the FASTA headers of the
// per-query <query name>.seqs.fa output files.  The header looks like this:
// >Score=614.319|Edges=1_2|Position=10|Alignment=369M|...
// Position is the 0-based offset of the alignment in the path traversal
// sequence; Alignment is a CIGAR-like string in amino acid space.
struct PathRacerAlignment {
    int position = 0;           // 0-based offset in the traversal
    int queryStart = 0;         // 1-based, in HMM (amino acid) coordinates
    int queryEnd = 0;           // 1-based, in HMM (amino acid) coordinates
    int alignedNtLength = 0;    // length of the aligned nucleotide span
};

} // namespace

// Parses the CIGAR-like alignment string.  Returns the number of leading and
// trailing deletions (these are HMM positions not aligned to anything) and
// the total nucleotide length of the aligned sequence span.
static bool parseAlignmentCigar(const QString &cigar,
                                int &leadingDeletions, int &trailingDeletions,
                                int &alignedNtLength) {
    static const QRegularExpression opRe("(\\d+)([MDI])");
    QRegularExpressionMatchIterator it = opRe.globalMatch(cigar);

    QVector<QPair<int, char>> ops;
    QString consumed;
    while (it.hasNext()) {
        auto match = it.next();
        ops.emplace_back(match.captured(1).toInt(), match.captured(2).at(0).toLatin1());
        consumed += match.captured(0);
    }

    // The whole string must be consumed by the regex
    if (ops.empty() || consumed != cigar)
        return false;

    leadingDeletions = 0;
    trailingDeletions = 0;
    alignedNtLength = 0;

    size_t i = 0;
    while (i < ops.size() && ops[i].second == 'D') {
        leadingDeletions += ops[i].first;
        ++i;
    }
    int j = int(ops.size()) - 1;
    while (j >= int(i) && ops[j].second == 'D') {
        trailingDeletions += ops[j].first;
        --j;
    }
    for (; i <= size_t(j); ++i) {
        // Both matches and insertions consume sequence (nucleotide) space
        if (ops[i].second == 'M' || ops[i].second == 'I')
            alignedNtLength += 3 * ops[i].first;
    }

    return true;
}

static std::vector<PathRacerPath>
parsePathRacerPaths(const QString &annotatedGfa) {
    std::vector<PathRacerPath> res;

    static const QRegularExpression pathLineRe("^P\\t(\\S+)\\t(\\S+)(?:\\t|$)");
    static const QRegularExpression pathNameRe("^(.*)_(\\d+)_score_(\\S+)$");

    for (const auto &line : annotatedGfa.split('\n', Qt::SkipEmptyParts)) {
        auto lineMatch = pathLineRe.match(line);
        if (!lineMatch.hasMatch())
            continue;

        auto nameMatch = pathNameRe.match(lineMatch.captured(1));
        if (!nameMatch.hasMatch())
            continue;

        bool idxOk = false, scoreOk = false;
        PathRacerPath path;
        path.queryName = nameMatch.captured(1);
        path.idx = nameMatch.captured(2).toInt(&idxOk);
        path.bitScore = nameMatch.captured(3).toDouble(&scoreOk);
        path.pathString = lineMatch.captured(2);
        if (idxOk && scoreOk)
            res.push_back(path);
    }

    return res;
}

// Parses all per-query <query name>.seqs.fa files from the PathRacer output
// directory.
static QMap<QString, QVector<PathRacerAlignment>>
parsePathRacerAlignments(const QDir &outputDir) {
    QMap<QString, QVector<PathRacerAlignment>> res;

    const auto seqsFiles =
            outputDir.entryInfoList(QStringList() << "*.seqs.fa", QDir::Files);
    for (const auto &fi : seqsFiles) {
        QString fileName = fi.fileName();
        if (!fileName.endsWith(".seqs.fa"))
            continue;
        QString queryName = fileName.left(fileName.length() - int(strlen(".seqs.fa")));

        QFile file(fi.absoluteFilePath());
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
            continue;

        QTextStream in(&file);
        while (!in.atEnd()) {
            QString header = in.readLine();
            if (!header.startsWith('>'))
                continue;

            // Parse the '|' separated key=value fields
            int position = -1, alignedNtLength = -1;
            int leadingDeletions = 0, trailingDeletions = 0;
            for (const auto &field : header.mid(1).split('|')) {
                int eq = field.indexOf('=');
                if (eq < 0)
                    continue;
                QString key = field.left(eq);
                QString value = field.mid(eq + 1);
                bool ok = false;
                if (key == "Position") {
                    position = value.toInt(&ok);
                    if (!ok) position = -1;
                } else if (key == "Alignment") {
                    if (!parseAlignmentCigar(value, leadingDeletions, trailingDeletions,
                                             alignedNtLength))
                        alignedNtLength = -1;
                }
            }

            if (position < 0 || alignedNtLength < 0)
                continue;

            PathRacerAlignment alignment;
            alignment.position = position;
            alignment.alignedNtLength = alignedNtLength;
            alignment.queryStart = leadingDeletions + 1;
            alignment.queryEnd = -trailingDeletions; // length is set below
            res[queryName].push_back(alignment);
        }
    }

    // Set the query ends using the HMM lengths
    // (done by the caller, which knows the queries)

    return res;
}

static std::pair<NodeHits, PathHits>
buildHitsFromPathRacerOutput(const QDir &outputDir, Queries &queries) {
    NodeHits nodeHits; PathHits pathHits;

    QFile gfaFile(outputDir.filePath("graph_with_hmm_paths.gfa"));
    if (!gfaFile.open(QIODevice::ReadOnly | QIODevice::Text))
        return {nodeHits, pathHits};

    QTextStream gfaIn(&gfaFile);
    const auto prPaths = parsePathRacerPaths(gfaIn.readAll());
    const auto alignments = parsePathRacerAlignments(outputDir);

    for (const auto &prPath : prPaths) {
        Query *query = queries.getQueryFromName(prPath.queryName);
        if (query == nullptr)
            continue;

        const auto &queryAlignments = alignments.value(prPath.queryName);
        if (prPath.idx < 0 || prPath.idx >= int(queryAlignments.size()))
            continue;

        PathRacerAlignment alignment = queryAlignments[prPath.idx];
        // queryEnd was initialized to -trailingDeletions; resolve it against
        // the HMM length: queryEnd = HMM length - trailing deletions.
        alignment.queryEnd = int(query->getLength()) + alignment.queryEnd;

        // Check the user-defined filters.  PathRacer does not report
        // e-values, so the e-value filter cannot be applied here.
        if (g_settings->bitScoreFilter.on &&
            prPath.bitScore < g_settings->bitScoreFilter)
            continue;

        int alignmentLength = alignment.queryEnd - alignment.queryStart + 1;
        if (g_settings->alignmentLengthFilter.on &&
            alignmentLength < g_settings->alignmentLengthFilter)
            continue;

        if (g_settings->queryCoverageFilter.on) {
            double hitCoveragePercentage = 100.0 * Hit::getQueryCoverageFraction(query,
                                                                                 alignment.queryStart, alignment.queryEnd);
            if (hitCoveragePercentage < g_settings->queryCoverageFilter)
                continue;
        }

        QString pathStringFailure;
        Path path = Path::makeFromString(prPath.pathString, *g_assemblyGraph,
                                         false, &pathStringFailure);
        if (path.isEmpty())
            continue;

        // Convert the 0-based offset in the traversal sequence into
        // 1-based path coordinates.
        int pathStart = alignment.position + 1;
        int pathEnd = alignment.position + alignment.alignedNtLength;

        for (const auto &entry : path.getNodeCovering(pathStart, pathEnd)) {
            DeBruijnNode *node = entry.first;
            int nodeStart = entry.second.mapped_range.from;
            int nodeEnd = entry.second.mapped_range.to;

            nodeHits.emplace_back(query,
                                  new Hit(query, node,
                                          -1, nodeEnd - nodeStart + 1,
                                          -1, -1,
                                          alignment.queryStart, alignment.queryEnd,
                                          nodeStart, nodeEnd,
                                          SciNot(0.0), prPath.bitScore));
        }
    }

    return { nodeHits, pathHits };
}

QString PathRacerSearch::doSearch(Queries &queries, QString extraParameters) {
    GraphSearchFinishedRAII watcher(this);

    m_lastError = "";
    if (!findTools())
        return m_lastError;

    // FIXME: Do we need proper mutex here?
    if (m_doSearch)
        return (m_lastError = "Search is already in progress");

    if (queries.empty())
        return (m_lastError = "There are no queries to search for");

    QTemporaryFile tmpQueryFile(temporaryDir().filePath("queries.XXXXXX.hmm"));
    if (!tmpQueryFile.open())
        return (m_lastError = "Failed to create temporary query file");

    writeQueryFile(&tmpQueryFile, queries);

    QString outputDirName = temporaryDir().filePath("pathracer_output");
    QDir outputDir(outputDirName);
    if (outputDir.exists())
        outputDir.removeRecursively();
    if (!QDir().mkpath(outputDirName))
        return (m_lastError = "Failed to create the PathRacer output directory");

    QStringList pathracerOptions;
    pathracerOptions << tmpQueryFile.fileName()
                     << temporaryDir().filePath("graph.gfa")
                     << "--output" << outputDirName
                     << extraParameters.split(" ", Qt::SkipEmptyParts);

    m_cancelSearch = false;
    m_doSearch = new QProcess();
    m_doSearch->start(m_pathracerCommand, pathracerOptions);

    bool finished = m_doSearch->waitForFinished(-1);

    if (m_doSearch->exitCode() != 0 || !finished) {
        if (m_cancelSearch) {
            m_lastError = "PathRacer search cancelled.";
        } else {
            m_lastError = "There was a problem running the PathRacer search";
            QString stdErr = m_doSearch->readAllStandardError();
            m_lastError += stdErr.isEmpty() ? "." : ":\n\n" + stdErr;
        }

        m_doSearch->deleteLater();
        m_doSearch = nullptr;

        return m_lastError;
    }

    m_doSearch->deleteLater();
    m_doSearch = nullptr;

    if (m_cancelSearch)
        return (m_lastError = "PathRacer search cancelled");

    auto [nodeHits, pathHits] = buildHitsFromPathRacerOutput(outputDir, queries);
    queries.addNodeHits(nodeHits);
    queries.findQueryPaths();
    queries.addPathHits(pathHits);
    queries.searchOccurred();

    m_lastError = "";

    return m_lastError;
}

QString PathRacerSearch::doAutoGraphSearch(const AssemblyGraph &graph, QString queriesFilename,
                                           bool includePaths,
                                           QString extraParameters) {
    cleanUp();

    QString maybeError = buildDatabase(graph, includePaths); // It is expected that buildDatabase will setup last error as well
    if (!maybeError.isEmpty())
        return maybeError;

    loadQueriesFromFile(queriesFilename);

    maybeError = doSearch(queries(), extraParameters);
    if (!maybeError.isEmpty())
        return maybeError;

    return "";
}

// This function returns the number of queries loaded from the HMM file.
int PathRacerSearch::loadQueriesFromFile(QString fullFileName) {
    m_lastError = "";

    int queriesBefore = int(getQueryCount());

    std::vector<QString> queryNames;
    std::vector<QByteArray> querySequences;
    std::vector<unsigned> queryLengths;
    std::vector<bool> queryProtHmms;
    if (!utils::readHmmFile(fullFileName,
                            queryNames, queryLengths, querySequences, queryProtHmms)) {
        m_lastError = "Failed to parse HMM file: " + fullFileName;
        return 0;
    }

    for (size_t i = 0; i < queryNames.size(); ++i)
        addQuery(new Query(cleanQueryName(queryNames[i]),
                           QString(queryLengths[i],
                                   queryProtHmms[i] ? 'F' : 'N'),
                           querySequences[i]));

    int queriesAfter = int(getQueryCount());
    return queriesAfter - queriesBefore;
}

void PathRacerSearch::cancelDatabaseBuild() {}

void PathRacerSearch::cancelSearch() {
    if (!m_doSearch)
        return;

    m_cancelSearch = true;
    if (m_doSearch)
        m_doSearch->kill();
}
