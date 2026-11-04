// Copyright 2022 Anton Korobeynikov

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

#include "graphscope.h"

#include "graphsearch/queries.h"

#include "program/settings.h"

#include "assemblygraph.h"
#include "path.h"

#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace graph {
    static QStringList splitCommaSeparated(QString text) {
        text = text.simplified();
        QStringList tokens = text.split(',');
        QStringList out;
        for (QString token : tokens) {
            token = token.trimmed();
            if (!token.isEmpty())
                out.push_back(token);
        }
        return out;
    }

    static void appendUnique(std::vector<DeBruijnNode *> &nodes,
                             std::unordered_set<DeBruijnNode *> &seen,
                             const std::vector<DeBruijnNode *> &more) {
        for (auto *node : more) {
            if (seen.insert(node).second)
                nodes.push_back(node);
        }
    }

    std::vector<DeBruijnNode *>
    getStartingNodes(QString *errorTitle, QString *errorMessage,
                     const AssemblyGraph &graph, const Scope &graphScope) {
        switch (graphScope.graphScope()) {
            default:
                return {};
            case AROUND_NODE: {
                if (AssemblyGraph::checkIfStringHasNodes(graphScope.nodeList())) {
                    *errorTitle = "No starting nodes";
                    *errorMessage = "Please enter at least one node when drawing the graph using the 'Around node(s)' scope. "
                                    "Separate multiple nodes with commas.";
                    return {};
                }

                // Make sure the nodes the user typed in are actually in the graph.
                std::vector<QString> nodesNotInGraph;
                std::vector<DeBruijnNode *> nodesInGraph =
                        graph.getNodesFromStringList(graphScope.nodeList(),
                                                     g_settings->startingNodesExactMatch,
                                                     &nodesNotInGraph);
                if (!nodesNotInGraph.empty()) {
                    *errorTitle = "Nodes not found";
                    *errorMessage =
                            AssemblyGraph::generateNodesNotFoundErrorMessage(nodesNotInGraph,
                                                                             g_settings->startingNodesExactMatch);
                }

                return nodesInGraph;
            }
            case AROUND_PATHS: {
                auto pathIt = graph.m_deBruijnGraphPaths.find(graphScope.path().toStdString());
                if (pathIt == graph.m_deBruijnGraphPaths.end()) {
                    *errorTitle = "Invalid path";
                    *errorMessage = "No path with such name is loaded";
                    return {};
                }

                return pathIt->nodes();
            }
            case AROUND_WALKS: {
                auto pathIt = graph.m_deBruijnGraphWalks.find(graphScope.walk().toStdString());
                if (pathIt == graph.m_deBruijnGraphWalks.end()) {
                    *errorTitle = "Invalid walk";
                    *errorMessage = "No walk with such sequence name is loaded";
                    return {};
                }

                return pathIt->walk.nodes();
            }
            case AROUND_COMPONENT: {
                QStringList pathNames = splitCommaSeparated(graphScope.componentPaths());
                QStringList walkNames = splitCommaSeparated(graphScope.componentWalks());
                bool hasNodes = !AssemblyGraph::checkIfStringHasNodes(graphScope.componentNodes());

                if (!hasNodes && pathNames.empty() && walkNames.empty()) {
                    *errorTitle = "No starting nodes";
                    *errorMessage = "Please enter at least one node, path or walk when drawing the graph using the "
                                    "'Around connected component' scope. Separate multiple names with commas.";
                    return {};
                }

                std::unordered_set<DeBruijnNode *> seen;
                std::vector<DeBruijnNode *> nodesInGraph;

                if (hasNodes) {
                    std::vector<QString> nodesNotInGraph;
                    auto fromNodes = graph.getNodesFromStringList(graphScope.componentNodes(),
                                                                  g_settings->startingNodesExactMatch,
                                                                  &nodesNotInGraph);
                    if (!nodesNotInGraph.empty()) {
                        *errorTitle = "Nodes not found";
                        *errorMessage =
                                AssemblyGraph::generateNodesNotFoundErrorMessage(nodesNotInGraph,
                                                                                 g_settings->startingNodesExactMatch);
                        return {};
                    }
                    appendUnique(nodesInGraph, seen, fromNodes);
                }

                QStringList missingPaths;
                for (const QString &name : pathNames) {
                    auto pathIt = graph.m_deBruijnGraphPaths.find(name.toStdString());
                    if (pathIt == graph.m_deBruijnGraphPaths.end()) {
                        missingPaths.push_back(name);
                        continue;
                    }
                    appendUnique(nodesInGraph, seen, pathIt->nodes());
                }
                if (!missingPaths.empty()) {
                    *errorTitle = "Invalid path";
                    *errorMessage = "No path with such name is loaded: " + missingPaths.join(", ");
                    return {};
                }

                QStringList missingWalks;
                for (const QString &name : walkNames) {
                    auto walkIt = graph.m_deBruijnGraphWalks.find(name.toStdString());
                    if (walkIt == graph.m_deBruijnGraphWalks.end()) {
                        missingWalks.push_back(name);
                        continue;
                    }
                    appendUnique(nodesInGraph, seen, walkIt->walk.nodes());
                }
                if (!missingWalks.empty()) {
                    *errorTitle = "Invalid walk";
                    *errorMessage = "No walk with such sequence name is loaded: " + missingWalks.join(", ");
                    return {};
                }

                if (nodesInGraph.empty()) {
                    *errorTitle = "No starting nodes";
                    *errorMessage = "The specified nodes, paths and walks did not match anything in the graph.";
                    return {};
                }

                return nodesInGraph;
            }
            case AROUND_QUERY_HITS: {
                std::vector<DeBruijnNode *> startingNodes;
                if (const auto *queries = graphScope.queries())
                    startingNodes = queries->getNodesFromHits(graphScope.queryName());
                if (startingNodes.empty()) {
                    *errorTitle = "No hits";
                    *errorMessage = "To draw the graph around graph hits, you must first conduct a graph search.";
                }
                return startingNodes;
            }
            case DEPTH_RANGE: {
                if (graphScope.minDepth() > graphScope.maxDepth()) {
                    *errorTitle = "Invalid depth range";
                    *errorMessage = "The maximum depth must be greater than or equal to the minimum depth.";
                    return {};
                }

                auto startingNodes = graph.getNodesInDepthRange(graphScope.minDepth(), graphScope.maxDepth());
                if (startingNodes.empty()) {
                    *errorTitle = "No nodes in range";
                    *errorMessage = "There are no nodes with depths in the specified range.";
                }
                return startingNodes;
            }
            case LENGTH_RANGE: {
                if (graphScope.minLength() > graphScope.maxLength()) {
                    *errorTitle = "Invalid length range";
                    *errorMessage = "The maximum length must be greater than or equal to the minimum length.";
                    return {};
                }

                auto startingNodes = graph.getNodesInLengthRange(graphScope.minLength(), graphScope.maxLength());
                if (startingNodes.empty()) {
                    *errorTitle = "No nodes in range";
                    *errorMessage = "There are no nodes with lengths in the specified range.";
                }
                return startingNodes;
            }
            case CIRCULAR_PATHS: {
                // GFA/SPAdes/GAF paths are always loaded as non-circular Path
                // objects (they do not include the closing edge), so we detect
                // circularity structurally: a path is circular if the graph
                // contains an edge leading from its last node back to its
                // first node.
                const bool includeSelfLoops = g_settings->circularPathsIncludeSelfLoops;
                // Set of pointers to all live nodes, used to detect paths
                // that still reference nodes deleted by graph edits (e.g.
                // node merging) without dereferencing them.
                std::unordered_set<const DeBruijnNode *> liveNodes;
                liveNodes.reserve(graph.m_deBruijnGraphNodes.size());
                for (const DeBruijnNode *node : graph.m_deBruijnGraphNodes)
                    liveNodes.insert(node);

                auto isCircularPath = [includeSelfLoops, &liveNodes](const Path &path) {
                    const auto &nodes = path.nodes();
                    if (nodes.empty())
                        return false;
                    for (const auto *node : nodes) {
                        if (liveNodes.find(node) == liveNodes.end())
                            return false;
                    }
                    // A single-node path is circular if the node connects
                    // back to itself (a self-loop edge).  Such paths are only
                    // considered when the 'include self-loop' setting is on.
                    if (nodes.size() == 1 && !includeSelfLoops)
                        return false;
                    const DeBruijnNode *lastNode = nodes.back();
                    const DeBruijnNode *firstNode = nodes.front();
                    for (const auto *edge : lastNode->edges()) {
                        if (edge->getStartingNode() == lastNode &&
                            edge->getEndingNode() == firstNode)
                            return true;
                    }
                    return false;
                };

                std::vector<DeBruijnNode *> startingNodes;
                std::unordered_set<DeBruijnNode *> seen;
                for (const auto &entry : graph.m_deBruijnGraphPaths) {
                    if (!isCircularPath(entry))
                        continue;
                    appendUnique(startingNodes, seen, entry.nodes());
                }

                if (startingNodes.empty()) {
                    *errorTitle = "No circular paths";
                    *errorMessage = "There are no circular paths in the graph. Note that circular paths are only "
                                    "detected when the last node of a path connects back to its first node.";
                }
                return startingNodes;
            }
            case CIRCULAR_CONTIGS: {
                // Unlike CIRCULAR_PATHS, this scope does not rely on path
                // records (which may be absent, e.g. in pre-scaffolding
                // SPAdes graphs).  Instead, circularity is detected purely
                // structurally: a node is considered circular if it has a
                // self-loop edge or if it lies on a directed cycle, i.e. it
                // belongs to a non-trivial strongly connected component.
                // An iterative version of Tarjan's SCC algorithm is used.

                auto hasSelfLoop = [](DeBruijnNode *node) {
                    for (const auto *edge : node->edges())
                        if (edge->getStartingNode() == node &&
                            edge->getEndingNode() == node)
                            return true;
                    return false;
                };

                auto successorsOf = [](DeBruijnNode *node) {
                    std::vector<DeBruijnEdge *> out;
                    for (auto *edge : node->edges())
                        if (edge->getStartingNode() == node)
                            out.push_back(edge);
                    return out;
                };

                std::unordered_map<DeBruijnNode *, int> indices, lowlinks;
                std::unordered_set<DeBruijnNode *> onStack;
                std::vector<DeBruijnNode *> sccStack;
                indices.reserve(graph.m_deBruijnGraphNodes.size());
                lowlinks.reserve(graph.m_deBruijnGraphNodes.size());
                int nextIndex = 0;

                std::unordered_set<DeBruijnNode *> cyclicNodes;

                struct Frame {
                    DeBruijnNode *node;
                    std::vector<DeBruijnEdge *> successors;
                    size_t next = 0;
                };
                std::vector<Frame> dfs;

                for (DeBruijnNode *root : graph.m_deBruijnGraphNodes) {
                    if (indices.count(root))
                        continue;
                    indices[root] = lowlinks[root] = nextIndex++;
                    sccStack.push_back(root);
                    onStack.insert(root);
                    dfs.push_back({root, successorsOf(root), 0});

                    while (!dfs.empty()) {
                        DeBruijnNode *node = dfs.back().node;
                        if (dfs.back().next < dfs.back().successors.size()) {
                            DeBruijnEdge *edge = dfs.back().successors[dfs.back().next++];
                            DeBruijnNode *succ = edge->getEndingNode();
                            if (!indices.count(succ)) {
                                indices[succ] = lowlinks[succ] = nextIndex++;
                                sccStack.push_back(succ);
                                onStack.insert(succ);
                                dfs.push_back({succ, successorsOf(succ), 0});
                            } else if (onStack.count(succ)) {
                                lowlinks[node] = std::min(lowlinks[node], indices[succ]);
                            }
                        } else {
                            dfs.pop_back();
                            if (!dfs.empty())
                                lowlinks[dfs.back().node] =
                                        std::min(lowlinks[dfs.back().node], lowlinks[node]);
                            if (lowlinks[node] == indices[node]) {
                                // Root of an SCC: pop the component off the stack.
                                std::vector<DeBruijnNode *> component;
                                while (true) {
                                    DeBruijnNode *w = sccStack.back();
                                    sccStack.pop_back();
                                    onStack.erase(w);
                                    component.push_back(w);
                                    if (w == node)
                                        break;
                                }
                                if (component.size() > 1)
                                    cyclicNodes.insert(component.begin(), component.end());
                                else if (hasSelfLoop(component[0]))
                                    cyclicNodes.insert(component[0]);
                            }
                        }
                    }
                }

                // Since edges come in reverse complement pairs, SCCs are
                // closed under reverse complementation.  Report only one
                // node of each {node, reverse complement} pair.
                std::vector<DeBruijnNode *> startingNodes;
                std::unordered_set<DeBruijnNode *> seen;
                for (DeBruijnNode *node : graph.m_deBruijnGraphNodes) {
                    if (!cyclicNodes.count(node))
                        continue;
                    DeBruijnNode *canonical = node->getCanonical();
                    if (!canonical)
                        canonical = node;
                    if (seen.insert(canonical).second)
                        startingNodes.push_back(canonical);
                }

                if (startingNodes.empty()) {
                    *errorTitle = "No circular contigs";
                    *errorMessage = "There are no circular contigs in the graph. Note that circular contigs are "
                                    "detected structurally: a node is circular if it has a self-loop edge or if "
                                    "it lies on a directed cycle in the graph.";
                }
                return startingNodes;
            }
        }

        return {};
    }

    Scope scope(GraphScope graphScope, const QString &nodesList,
                double minDepthRange, double maxDepthRange,
                const search::Queries *queries, const QString &queryName,
                const QString &pathName, unsigned distance,
                const QString &walkName,
                long long minLengthRange, long long maxLengthRange) {
        switch (graphScope) {
            case WHOLE_GRAPH:
                return Scope::wholeGraph();
            case AROUND_NODE:
                return Scope::aroundNodes(nodesList, distance);
            case AROUND_PATHS:
                return Scope::aroundPath(pathName, distance);
            case AROUND_WALKS:
                return Scope::aroundWalk(walkName.isEmpty() ? pathName : walkName, distance);
            case AROUND_QUERY_HITS:
                return Scope::aroundHits(queries, queryName, distance);
            case DEPTH_RANGE:
                return Scope::depthRange(minDepthRange, maxDepthRange);
            case LENGTH_RANGE:
                return Scope::lengthRange(minLengthRange, maxLengthRange);
            case CIRCULAR_PATHS:
                return Scope::circularPaths();
            case CIRCULAR_CONTIGS:
                return Scope::circularContigs();
            case AROUND_COMPONENT:
                return Scope::aroundComponent(nodesList, pathName, walkName);
        }

        assert(0 && "Invalid scope!");
        return Scope::wholeGraph();
    }

    Scope Scope::aroundHits(const search::Queries &queries, const QString& queryName,
                            unsigned distance) {
        return aroundHits(&queries, queryName, distance);
    }

    Scope Scope::aroundHits(const search::Queries *queries, const QString& queryName,
                            unsigned distance) {
        Scope res;
        res.m_scope = AROUND_QUERY_HITS;
        res.m_opt = std::make_pair(queries, queryName);
        res.m_distance = distance;

        return res;
    }

    Scope Scope::aroundComponent(QString nodes, QString paths, QString walks) {
        Scope res;
        res.m_scope = AROUND_COMPONENT;
        res.m_opt = ComponentSeed{std::move(nodes), std::move(paths), std::move(walks)};
        return res;
    }
}
