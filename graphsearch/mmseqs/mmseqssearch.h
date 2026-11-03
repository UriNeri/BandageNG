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

#pragma once

#include "graphsearch/graphsearch.h"

#include <QDir>
#include <QString>

class QProcess;

namespace search {

class Queries;

// Searcher based on MMseqs2 (https://github.com/soedinglab/MMseqs2).
// It uses the "easy-search" workflow which accepts plain FASTA files
// and outputs BLAST tabular format, so no explicit database build
// step is necessary.
class MmseqsSearch : public search::GraphSearch {
    Q_OBJECT
public:
    explicit MmseqsSearch(const QDir &workDir = QDir::temp(), QObject *parent = nullptr);
    virtual ~MmseqsSearch() = default;

    QString doAutoGraphSearch(const AssemblyGraph &graph, QString queriesFilename,
                              bool includePaths = false,
                              QString extraParameters = "") override;
    int loadQueriesFromFile(QString fullFileName) override;
    QString buildDatabase(const AssemblyGraph &graph,
                          bool includePaths = true) override;
    QString doSearch(QString extraParameters) override;
    QString doSearch(search::Queries &queries, QString extraParameters) override;

    QString name() const override { return "MMseqs2"; }
    QString queryFormat() const override { return "FASTA"; }
    QString annotationGroupName() const override { return "MMseqs2 hits"; };

public slots:
    void cancelDatabaseBuild() override;
    void cancelSearch() override;

private:
    bool findTools();

    bool m_cancelBuildDatabase = false, m_cancelSearch = false;

    QProcess *m_doSearch = nullptr;
    QString m_mmseqsCommand;
};

}
