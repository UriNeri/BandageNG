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

#pragma once

#include "graphsearch/graphsearch.h"

#include <QDir>
#include <QString>

class QProcess;

namespace search {

class Queries;

// Graph search backend based on PathRacer
// (https://github.com/ablab/spades/tree/pathracer), the SPAdes tool for
// HMM-to-graph alignment.  PathRacer aligns protein (or nucleotide) HMMs
// directly to the assembly graph, performing its own frame-shift aware
// six-frame translation internally, so Bandage does not need to do any
// translation itself.
class PathRacerSearch : public GraphSearch {
    Q_OBJECT
public:
    explicit PathRacerSearch(const QDir &workDir = QDir::temp(), QObject *parent = nullptr);
    virtual ~PathRacerSearch() = default;

    QString doAutoGraphSearch(const AssemblyGraph &graph, QString queriesFilename,
                              bool includePaths = false,
                              QString extraParameters = "") override;
    int loadQueriesFromFile(QString fullFileName) override;
    QString buildDatabase(const AssemblyGraph &graph,
                          bool includePaths = true) override;
    QString doSearch(QString extraParameters) override;
    QString doSearch(search::Queries &queries, QString extraParameters) override;

    bool allowManualQueries() const override { return false; }
    QString name() const override { return "PathRacer"; }
    QString queryFormat() const override { return "HMM"; }
    QString annotationGroupName() const override { return "PathRacer hits"; }

public slots:
    void cancelDatabaseBuild() override;
    void cancelSearch() override;

private:
    bool findTools();

    bool m_cancelSearch = false;

    QProcess *m_doSearch = nullptr;
    QString m_pathracerCommand;
};

}
