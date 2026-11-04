// Copyright 2023 Anton Korobeynikov

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


#include "diamondsearch.h"

#include "graphsearch/graphsearch.h"
#include "program/globals.h"
#include "program/settings.h"

#include "graph/assemblygraph.h"
#include "graph/sequenceutils.h"
#include "io/fileutils.h"

#include <QDir>
#include <QProcess>
#include <QTemporaryFile>

#include <array>
#include <map>

using namespace search;

DiamondSearch::DiamondSearch(const QDir &workDir, QObject *parent)
        : GraphSearch(workDir, parent) {}

bool DiamondSearch::findTools() {
    if (!findProgram("diamond", g_settings->diamondCommand, &m_diamondCommand, "DIAMOND"))
        return false;

    return true;
}

namespace {

// Standard genetic code translation table (codon -> amino acid).
// 'X' is used for stop codons and ambiguous codons.
char translateCodon(const char *codon) {
    static const std::map<std::string, char> table = {
        {"TTT", 'F'}, {"TTC", 'F'}, {"TTA", 'L'}, {"TTG", 'L'},
        {"CTT", 'L'}, {"CTC", 'L'}, {"CTA", 'L'}, {"CTG", 'L'},
        {"ATT", 'I'}, {"ATC", 'I'}, {"ATA", 'I'}, {"ATG", 'M'},
        {"GTT", 'V'}, {"GTC", 'V'}, {"GTA", 'V'}, {"GTG", 'V'},
        {"TCT", 'S'}, {"TCC", 'S'}, {"TCA", 'S'}, {"TCG", 'S'},
        {"CCT", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'},
        {"ACT", 'T'}, {"ACC", 'T'}, {"ACA", 'T'}, {"ACG", 'T'},
        {"GCT", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'},
        {"TAT", 'Y'}, {"TAC", 'Y'}, {"TAA", 'X'}, {"TAG", 'X'},
        {"CAT", 'H'}, {"CAC", 'H'}, {"CAA", 'Q'}, {"CAG", 'Q'},
        {"AAT", 'N'}, {"AAC", 'N'}, {"AAA", 'K'}, {"AAG", 'K'},
        {"GAT", 'D'}, {"GAC", 'D'}, {"GAA", 'E'}, {"GAG", 'E'},
        {"TGT", 'C'}, {"TGC", 'C'}, {"TGA", 'X'}, {"TGG", 'W'},
        {"CGT", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'},
        {"AGT", 'S'}, {"AGC", 'S'}, {"AGA", 'R'}, {"AGG", 'R'},
        {"GGT", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'},
    };

    auto it = table.find(std::string(codon, 3));
    return it == table.end() ? 'X' : it->second;
}

char complementBase(char c) {
    switch (c) {
        case 'A': case 'a': return 'T';
        case 'T': case 't': return 'A';
        case 'G': case 'g': return 'C';
        case 'C': case 'c': return 'G';
        default: return 'N';
    }
}

QByteArray reverseComplement(const QByteArray &seq) {
    QByteArray res;
    res.reserve(seq.size());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it)
        res.append(complementBase(*it));
    return res;
}

// Translate a nucleotide sequence in a single frame.
// frame: 0, 1, 2 for the forward strand offsets; -1, -2, -3 for the
// reverse strand offsets (i.e. translation starts at |frame|-1 bases
// from the end of the reverse complement).
QByteArray translateFrame(const QByteArray &seq, int frame) {
    QByteArray working;
    if (frame > 0) {
        working = seq.mid(frame - 1);
    } else {
        working = reverseComplement(seq).mid(-frame - 1);
    }

    QByteArray protein;
    protein.reserve(working.size() / 3);
    for (int i = 0; i + 2 < working.size(); i += 3)
        protein.append(translateCodon(working.constData() + i));

    return protein;
}

// The subject names used in the DIAMOND database encode the original
// sequence name and the frame: "<name>#f<frame>".
QString frameSubjectName(const QString &name, int frame) {
    return name + "#f" + QString::number(frame);
}

struct SubjectLocation {
    QString name;
    int frame = 0;
};

SubjectLocation parseSubjectName(const QString &subject) {
    int idx = subject.lastIndexOf("#f");
    if (idx < 0)
        return {subject, 0};

    bool ok = false;
    int frame = subject.mid(idx + 2).toInt(&ok);
    if (!ok)
        return {subject, 0};

    return {subject.left(idx), frame};
}

// Convert a subject (amino acid) coordinate to a nucleotide coordinate
// within the original sequence, taking the frame into account.
// The returned coordinate is 1-based.
int subjectPosToNuclPos(int aaPos, int frame, int seqLength) {
    // aaPos is 1-based within the translated frame.
    // For forward frames, nucleotide position = 3 * (aaPos - 1) + frame.
    // For reverse frames, count from the end of the sequence.
    int offset = 3 * (aaPos - 1);
    if (frame > 0)
        return frame + offset;
    else
        return seqLength - (-frame - 1) - offset;
}

} // namespace

// Build a protein database from the 6-frame translation of all graph
// sequences.  DIAMOND's makedb works on FASTA input.
QString DiamondSearch::buildDatabase(const AssemblyGraph &graph, bool includePaths) {
    DbBuildFinishedRAII watcher(this);

    m_lastError = "";
    if (!findTools())
        return m_lastError;

    if (m_buildDb)
        return (m_lastError = "Building is already in progress");

    m_cancelBuildDatabase = false;

    QFile file(temporaryDir().filePath("all_nodes_translated.fasta"));
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        return (m_lastError = "Failed to open: " + file.fileName());

    {
        QTextStream out(&file);
        auto writeTranslated = [&](const QString &name, const QByteArray &sequence) {
            if (m_cancelBuildDatabase)
                return false;

            if (sequence.isEmpty())
                return true;

            for (int frame : {1, 2, 3, -1, -2, -3}) {
                QByteArray protein = translateFrame(sequence, frame);
                if (protein.isEmpty())
                    continue;

                out << '>' << frameSubjectName(name, frame) << '\n'
                    << protein << '\n';
            }

            return true;
        };

        for (const auto *node : graph.m_deBruijnGraphNodes) {
            if (node->sequenceIsMissing())
                continue;

            // Only translate positive nodes - negative nodes are their
            // reverse complements and would duplicate the frames.
            if (!node->isPositiveNode())
                continue;

            QByteArray sequence = utils::sequenceToQByteArray(node->getSequence());
            sequence = sequence.simplified();
            sequence.replace(' ', "");
            if (!writeTranslated(node->getName(), sequence))
                return (m_lastError = "Build cancelled.");
        }

        if (includePaths) {
            for (auto it = graph.m_deBruijnGraphPaths.begin(); it != graph.m_deBruijnGraphPaths.end(); ++it) {
                QByteArray sequence = it.value().getPathSequence();
                if (!writeTranslated(it.key().c_str(), sequence))
                    return (m_lastError = "Build cancelled.");
            }
        }
    }

    // Make sure the graph has sequences
    bool atLeastOneSequence = false;
    for (const auto *node : graph.m_deBruijnGraphNodes) {
        if (!node->sequenceIsMissing()) {
            atLeastOneSequence = true;
            break;
        }
    }

    if (!atLeastOneSequence)
        return (m_lastError = "Cannot build the DIAMOND database as this graph contains no sequences");

    QStringList makeDbOptions;
    makeDbOptions << "makedb"
                  << "--in" << temporaryDir().filePath("all_nodes_translated.fasta")
                  << "--db" << temporaryDir().filePath("all_nodes.dmnd");

    m_buildDb = new QProcess();
    m_buildDb->start(m_diamondCommand, makeDbOptions);

    bool finished = m_buildDb->waitForFinished(-1);
    if (m_buildDb->exitCode() != 0 || !finished) {
        m_lastError = "There was a problem building the DIAMOND database";
        QString stdErr = m_buildDb->readAllStandardError();
        m_lastError += stdErr.isEmpty() ? "." : ":\n\n" + stdErr;
    } else if (m_cancelBuildDatabase)
        m_lastError = "Build cancelled.";
    else
        m_lastError = "";

    m_buildDb->deleteLater();
    m_buildDb = nullptr;

    return m_lastError;
}

QString DiamondSearch::doSearch(QString extraParameters) {
    return doSearch(queries(), extraParameters);
}

static void writeQueryFile(QFile *file,
                           const Queries &queries,
                           QuerySequenceType sequenceType) {
    QTextStream out(file);
    for (const auto *query: queries.queries()) {
        if (query->getSequenceType() != sequenceType)
            continue;

        out << '>' << query->getName() << '\n'
            << query->getSequence()
            << '\n';
    }
}

QString DiamondSearch::runOneDiamondSearch(QuerySequenceType sequenceType,
                                           const Queries &queries,
                                           const QString &extraParameters,
                                           bool &success) {
    QTemporaryFile tmpFile(temporaryDir().filePath(sequenceType == NUCLEOTIDE ?
                                                   "nucl_queries.XXXXXX.fasta" : "prot_queries.XXXXXX.fasta"));
    if (!tmpFile.open()) {
        m_lastError = "Failed to create temporary query file";
        success = false;
        return "";
    }

    writeQueryFile(&tmpFile, queries, sequenceType);

    // Nucleotide queries are searched against the translated database
    // with blastx; protein queries with blastp.
    QStringList diamondOptions;
    diamondOptions << (sequenceType == NUCLEOTIDE ? "blastx" : "blastp")
                   << "--query" << tmpFile.fileName()
                   << "--db" << temporaryDir().filePath("all_nodes.dmnd")
                   << "--outfmt" << "6";
    diamondOptions << extraParameters.split(" ", Qt::SkipEmptyParts);

    m_doSearch = new QProcess();
    m_doSearch->start(m_diamondCommand, diamondOptions);

    bool finished = m_doSearch->waitForFinished(-1);
    if (m_doSearch->exitCode() != 0 || !finished) {
        if (m_cancelSearch) {
            m_lastError = "DIAMOND search cancelled.";
        } else {
            m_lastError = "There was a problem running the DIAMOND search";
            QString stdErr = m_doSearch->readAllStandardError();
            m_lastError += stdErr.isEmpty() ? "." : ":\n\n" + stdErr;
        }
        m_doSearch->deleteLater();
        m_doSearch = nullptr;

        success = false;
        return m_lastError;
    }

    QString diamondOutput = m_doSearch->readAllStandardOutput();
    m_doSearch->deleteLater();
    m_doSearch = nullptr;

    success = true;
    return diamondOutput;
}

// DIAMOND outputs BLAST tabular format.  Subject coordinates are in
// amino acids and need to be converted back to nucleotide coordinates
// using the frame encoded in the subject name.
static std::pair<NodeHits, PathHits>
buildHitsFromDiamondOutput(QString diamondOutput,
                           Queries &queries);
static std::pair<NodeHits, PathHits>
buildHitsFromDiamondOutput(QString diamondOutput,
                           Queries &queries) {
    NodeHits nodeHits; PathHits pathHits;

    QStringList diamondHitList = diamondOutput.split("\n", Qt::SkipEmptyParts);

    for (const auto &hitString : diamondHitList) {
        QStringList alignmentParts = hitString.split('\t');

        if (alignmentParts.size() < 12)
            continue;

        QString queryName = alignmentParts[0];
        QString subjectLabel = alignmentParts[1];
        double percentIdentity = alignmentParts[2].toDouble();
        int alignmentLength = alignmentParts[3].toInt();
        int numberMismatches = alignmentParts[4].toInt();
        int numberGapOpens = alignmentParts[5].toInt();
        int queryStart = alignmentParts[6].toInt();
        int queryEnd = alignmentParts[7].toInt();
        int subjectStart = alignmentParts[8].toInt();
        int subjectEnd = alignmentParts[9].toInt();
        SciNot eValue(alignmentParts[10]);
        double bitScore = alignmentParts[11].toDouble();

        Query *query = queries.getQueryFromName(queryName);
        if (query == nullptr)
            continue;

        // Check the user-defined filters.
        if (g_settings->identityFilter.on &&
            percentIdentity < g_settings->identityFilter)
            continue;

        if (g_settings->evalueFilter.on &&
            eValue > g_settings->evalueFilter)
            continue;

        if (g_settings->bitScoreFilter.on &&
            bitScore < g_settings->bitScoreFilter)
            continue;

        if (g_settings->alignmentLengthFilter.on &&
            alignmentLength < g_settings->alignmentLengthFilter)
            continue;

        if (g_settings->queryCoverageFilter.on) {
            double hitCoveragePercentage = 100.0 * Hit::getQueryCoverageFraction(query,
                                                                                 queryStart, queryEnd);
            if (hitCoveragePercentage < g_settings->queryCoverageFilter)
                continue;
        }

        auto [nodeLabel, frame] = parseSubjectName(subjectLabel);

        // Convert the amino acid subject coordinates to nucleotide
        // coordinates.  We need the sequence length for reverse frames,
        // so look the node/path up first.
        auto nodeIt = g_assemblyGraph->m_deBruijnGraphNodes.find(nodeLabel.toStdString());
        auto pathIt = g_assemblyGraph->m_deBruijnGraphPaths.find(nodeLabel.toStdString());

        int seqLength = 0;
        if (nodeIt != g_assemblyGraph->m_deBruijnGraphNodes.end())
            seqLength = nodeIt.value()->getLength();
        else if (pathIt != g_assemblyGraph->m_deBruijnGraphPaths.end())
            seqLength = pathIt.value().getLength();
        else
            continue;

        int nodeStart = subjectPosToNuclPos(std::min(subjectStart, subjectEnd), frame, seqLength);
        int nodeEnd = subjectPosToNuclPos(std::max(subjectStart, subjectEnd), frame, seqLength);

        // Only save hits that are on forward strands.
        if (frame < 0)
            continue;

        if (nodeIt != g_assemblyGraph->m_deBruijnGraphNodes.end()) {
            nodeHits.emplace_back(query,
                                  new Hit(query, nodeIt.value(),
                                          percentIdentity, alignmentLength * 3,
                                          numberMismatches * 3, numberGapOpens,
                                          queryStart, queryEnd,
                                          nodeStart, nodeEnd, eValue, bitScore));
        }

        if (pathIt != g_assemblyGraph->m_deBruijnGraphPaths.end()) {
            pathHits.emplace_back(query, &pathIt.value(),
                                  Path::MappingRange{queryStart, queryEnd,
                                                     nodeStart, nodeEnd});
        }

    }

    return { nodeHits, pathHits };
}

QString DiamondSearch::doSearch(Queries &queries, QString extraParameters) {
    GraphSearchFinishedRAII watcher(this);

    m_lastError = "";
    if (!findTools())
        return m_lastError;

    // FIXME: Do we need proper mutex here?
    if (m_doSearch)
        return (m_lastError = "Search is already in progress");

    QString diamondOutput;
    bool success = false;
    if (queries.getQueryCount(NUCLEOTIDE) > 0 && !m_cancelSearch) {
        diamondOutput += runOneDiamondSearch(NUCLEOTIDE, queries, extraParameters, success);
        if (!success)
            return m_lastError;
    }

    if (queries.getQueryCount(PROTEIN) > 0 && !m_cancelSearch) {
        diamondOutput += runOneDiamondSearch(PROTEIN, queries, extraParameters, success);
        if (!success)
            return m_lastError;
    }

    if (m_cancelSearch)
        return (m_lastError = "DIAMOND search cancelled");

    auto [nodeHits, pathHits] = buildHitsFromDiamondOutput(diamondOutput, queries);
    queries.addNodeHits(nodeHits);
    queries.findQueryPaths();
    queries.addPathHits(pathHits);
    queries.searchOccurred();

    m_lastError = "";

    return m_lastError;
}

QString DiamondSearch::doAutoGraphSearch(const AssemblyGraph &graph, QString queriesFilename,
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

//This function returns the number of queries loaded from the FASTA file.
int DiamondSearch::loadQueriesFromFile(QString fullFileName) {
    m_lastError = "";
    int queriesBefore = int(getQueryCount());

    std::vector<QString> queryNames;
    std::vector<QByteArray> querySequences;
    if (!utils::readFastxFile(fullFileName, queryNames, querySequences)) {
        m_lastError = "Failed to parse FASTA file: " + fullFileName;
        return 0;
    }

    for (size_t i = 0; i < queryNames.size(); ++i) {
        //We only use the part of the query name up to the first space.
        QStringList queryNameParts = queryNames[i].split(" ");
        QString queryName;
        if (!queryNameParts.empty())
            queryName = cleanQueryName(queryNameParts[0]);

        addQuery(new Query(queryName, querySequences[i]));
    }

    int queriesAfter = int(getQueryCount());
    return queriesAfter - queriesBefore;
}

void DiamondSearch::cancelDatabaseBuild() {
    if (!m_buildDb)
        return;

    m_cancelBuildDatabase = true;
    if (m_buildDb)
        m_buildDb->kill();
}

void DiamondSearch::cancelSearch() {
    if (!m_doSearch)
        return;

    m_cancelSearch = true;
    if (m_doSearch)
        m_doSearch->kill();
}
