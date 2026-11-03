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


#include "mmseqssearch.h"

#include "graphsearch/graphsearch.h"
#include "program/settings.h"

#include "graph/assemblygraph.h"
#include "io/fileutils.h"

#include <QDir>
#include <QProcess>
#include <QTemporaryFile>

using namespace search;

MmseqsSearch::MmseqsSearch(const QDir &workDir, QObject *parent)
        : GraphSearch(workDir, parent) {}

bool MmseqsSearch::findTools() {
    if (!findProgram("mmseqs", &m_mmseqsCommand)) {
        m_lastError = "Error: The program mmseqs was not found.  Please install MMseqs2 to use this feature.";
        return false;
    }

    return true;
}

// MMseqs2 easy-search works directly on FASTA files, so "building the
// database" simply means writing all node (and optionally path) sequences
// into a single FASTA file.
QString MmseqsSearch::buildDatabase(const AssemblyGraph &graph, bool includePaths) {
    DbBuildFinishedRAII watcher(this);
    m_lastError = "";
    if (!findTools())
        return m_lastError;

    m_cancelBuildDatabase = false;

    QFile file(temporaryDir().filePath("all_nodes.fasta"));
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        return (m_lastError = "Failed to open: " + file.fileName());

    {
        QTextStream out(&file);
        for (const auto *node : graph.m_deBruijnGraphNodes) {
            if (m_cancelBuildDatabase)
                return (m_lastError = "Build cancelled.");

            out << node->getFasta(true, false, false);
        }

        if (includePaths) {
            for (auto it = graph.m_deBruijnGraphPaths.begin(); it != graph.m_deBruijnGraphPaths.end(); ++it) {
                if (m_cancelBuildDatabase)
                    return (m_lastError = "Build cancelled.");

                out << it.value().getFasta(it.key().c_str());
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
        return (m_lastError = "Cannot build the MMseqs2 database as this graph contains no sequences");

    return m_lastError;
}

QString MmseqsSearch::doSearch(QString extraParameters) {
    return doSearch(queries(), extraParameters);
}

static void writeQueryFile(QFile *file,
                           const Queries &queries) {
    QTextStream out(file);
    for (const auto *query: queries.queries()) {
        out << '>' << query->getName() << '\n'
            << query->getSequence()
            << '\n';
    }
}

// FIXME: Factor out common functionality
static QString getNodeNameFromString(const QString &nodeString) {
    QStringList nodeStringParts = nodeString.split("_");

    // The node string format should look like this:
    // NODE_nodename_length_123_cov_1.23
    if (nodeStringParts.size() < 6)
        return "";

    if (nodeStringParts.size() == 6)
        return nodeStringParts[1];

    // If the code got here, there are more than 6 parts.  This means there are
    // underscores in the node name (happens a lot with Trinity graphs).  So we
    // need to pull out the parts which constitute the name.
    int underscoreCount = nodeStringParts.size() - 6;
    QString nodeName = "";
    for (int i = 0; i <= underscoreCount; ++i) {
        nodeName += nodeStringParts[1+i];
        if (i < underscoreCount)
            nodeName += "_";
    }

    return nodeName;
}

// MMseqs2 easy-search outputs BLAST tabular format (outfmt 6 compatible):
// query, target, pident, alnlen, mismatch, gapopen, qstart, qend,
// tstart, tend, evalue, bits
static std::pair<NodeHits, PathHits>
buildHitsFromMmseqsOutput(QString mmseqsOutput,
                          Queries &queries) {
    NodeHits nodeHits; PathHits pathHits;

    QStringList mmseqsHitList = mmseqsOutput.split("\n", Qt::SkipEmptyParts);

    for (const auto &hitString : mmseqsHitList) {
        QStringList alignmentParts = hitString.split('\t');

        if (alignmentParts.size() < 12)
            continue;

        QString queryName = alignmentParts[0];
        QString nodeLabel = alignmentParts[1];
        double percentIdentity = alignmentParts[2].toDouble();
        int alignmentLength = alignmentParts[3].toInt();
        int numberMismatches = alignmentParts[4].toInt();
        int numberGapOpens = alignmentParts[5].toInt();
        int queryStart = alignmentParts[6].toInt();
        int queryEnd = alignmentParts[7].toInt();
        int nodeStart = alignmentParts[8].toInt();
        int nodeEnd = alignmentParts[9].toInt();
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

        auto nodeIt = g_assemblyGraph->m_deBruijnGraphNodes.find(getNodeNameFromString(nodeLabel).toStdString());
        if (nodeIt != g_assemblyGraph->m_deBruijnGraphNodes.end()) {
            // Only save hits that are on forward strands.
            if (nodeStart > nodeEnd)
                continue;

            nodeHits.emplace_back(query,
                                  new Hit(query, nodeIt.value(),
                                          percentIdentity, alignmentLength,
                                          numberMismatches, numberGapOpens,
                                          queryStart, queryEnd,
                                          nodeStart, nodeEnd, eValue, bitScore));
        }

        auto pathIt = g_assemblyGraph->m_deBruijnGraphPaths.find(nodeLabel.toStdString());
        if (pathIt != g_assemblyGraph->m_deBruijnGraphPaths.end()) {
            pathHits.emplace_back(query, &pathIt.value(),
                                  Path::MappingRange{queryStart, queryEnd,
                                                     nodeStart, nodeEnd});
        }

    }

    return { nodeHits, pathHits };
}

QString MmseqsSearch::doSearch(Queries &queries, QString extraParameters) {
    GraphSearchFinishedRAII watcher(this);

    m_lastError = "";
    if (!findTools())
        return m_lastError;

    // FIXME: Do we need proper mutex here?
    if (m_doSearch)
        return (m_lastError = "Search is already in progress");

    QTemporaryFile tmpFile(temporaryDir().filePath("queries.XXXXXX.fasta"));
    if (!tmpFile.open())
        return (m_lastError = "Failed to create temporary query file");

    writeQueryFile(&tmpFile, queries);

    // MMseqs2 easy-search: query, target, output file, temp dir.
    // The output goes to stdout in BLAST tabular format.
    QStringList mmseqsOptions;
    mmseqsOptions << "easy-search"
                  << tmpFile.fileName()
                  << temporaryDir().filePath("all_nodes.fasta")
                  << "stdout"
                  << temporaryDir().filePath("tmp")
                  << "--format-output" << "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits";
    mmseqsOptions << extraParameters.split(" ", Qt::SkipEmptyParts);

    m_cancelSearch = false;
    m_doSearch = new QProcess();
    m_doSearch->start(m_mmseqsCommand, mmseqsOptions);

    bool finished = m_doSearch->waitForFinished(-1);

    if (m_doSearch->exitCode() != 0 || !finished) {
        if (m_cancelSearch) {
            m_lastError = "MMseqs2 search cancelled.";
        } else {
            m_lastError = "There was a problem running the MMseqs2 search";
            QString stdErr = m_doSearch->readAllStandardError();
            m_lastError += stdErr.isEmpty() ? "." : ":\n\n" + stdErr;
        }

        m_doSearch->deleteLater();
        m_doSearch = nullptr;

        return m_lastError;
    }

    QString mmseqsOutput = m_doSearch->readAllStandardOutput();
    m_doSearch->deleteLater();
    m_doSearch = nullptr;

    if (m_cancelSearch)
        return (m_lastError = "MMseqs2 search cancelled");

    auto [nodeHits, pathHits] = buildHitsFromMmseqsOutput(mmseqsOutput, queries);
    queries.addNodeHits(nodeHits);
    queries.findQueryPaths();
    queries.addPathHits(pathHits);
    queries.searchOccurred();

    m_lastError = "";

    return m_lastError;
}

QString MmseqsSearch::doAutoGraphSearch(const AssemblyGraph &graph, QString queriesFilename,
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
int MmseqsSearch::loadQueriesFromFile(QString fullFileName) {
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

void MmseqsSearch::cancelDatabaseBuild() {
    m_cancelBuildDatabase = true;
}

void MmseqsSearch::cancelSearch() {
    if (!m_doSearch)
        return;

    m_cancelSearch = true;
    if (m_doSearch)
        m_doSearch->kill();
}
