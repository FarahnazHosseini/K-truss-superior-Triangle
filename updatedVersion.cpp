#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <time.h>
using namespace std;

class graph {
public:
    int totalEdge;
    int totalNode;
    vector<int> csr;
    vector< vector <int > > supports;
    vector<int> degree;
    vector <pair <int, int> > sourceEdgeSet;
    vector < vector <int> > superiorTriangles;
    vector <vector <int> > neighborList;
    vector < vector <int> > results;
    vector < pair <int, int> > uniqueEdges;

    void generateneighborList();
    void findUnionEdges();
    void addEdgeToNeighborList(int, int);
    void findKinUnionEdges();
    void setGraph(int, int);
    void optimizedKtrussForUnionEdges();
    void Ktruss();
    void KtrussNeighbourList();
    void findSuperiorTriangles();
    void printSupport();
    void findCommonNeighbors();
    void printGraphVector();
    vector<int> getNeighbors(int node);
    void processTriangleEdges(const vector<int>& triangle);
    void findSourceEdgeSet(int a, int b);
    friend void readFile(graph& g, string fileName);
    friend void printGraph(graph g);
    friend void printEdge(graph g);

};
//...........................
vector <int> graph::getNeighbors(int node) {
    vector <int> neighbors;
    int start = csr[node];
    int end = (node == totalNode - 1) ? csr.size() : csr[node + 1];
    for (int i = start; i < end; ++i) {
        neighbors.push_back(csr[i]);
    }
    return neighbors;
}
//............................
void addEdge(vector<pair <int, int> >& sourceEdgeSet, int u, int v) {
    pair<int, int> edge = make_pair(min(u, v), max(u, v));
    for (size_t i = 0; i < sourceEdgeSet.size(); ++i) {
        if (sourceEdgeSet[i] == edge) {
            return; // Avoid adding duplicates
        }
    }
    sourceEdgeSet.push_back(edge);
    ////add new edge to neighbor list

}
//...........................
// Helper function to check if an edge exists in the vector
bool edgeExists(const vector <pair <int, int> >& edgeList, const pair <int, int>& edge) {
    for (int i = 0; i < edgeList.size(); i++) {
        if (edgeList[i] == edge) {
            return true;
        }
    }
    return false;
}
//............................................
// Function to find the union of edges in superior triangles and source edge set
void graph::findUnionEdges() {
    // Add edges from superior triangles
    for (int i = 0; i < superiorTriangles.size(); i++) {
        vector<int> triangle = superiorTriangles[i];
        pair<int, int> edge1 = make_pair(min(triangle[0], triangle[1]), max(triangle[0], triangle[1]));
        pair<int, int> edge2 = make_pair(min(triangle[1], triangle[2]), max(triangle[1], triangle[2]));
        pair<int, int> edge3 = make_pair(min(triangle[0], triangle[2]), max(triangle[0], triangle[2]));
        if (!edgeExists(uniqueEdges, edge1)) {
            uniqueEdges.push_back(edge1);
        }
        if (!edgeExists(uniqueEdges, edge2)) {
            uniqueEdges.push_back(edge2);
        }
        if (!edgeExists(uniqueEdges, edge3)) {
            uniqueEdges.push_back(edge3);
        }
    }

    // Add edges from the source edge set
    for (int i = 0; i < sourceEdgeSet.size(); i++) {
        pair<int, int> edge = sourceEdgeSet[i];
        pair<int, int> sortedEdge = make_pair(min(edge.first, edge.second), max(edge.first, edge.second));

        if (!edgeExists(uniqueEdges, sortedEdge)) {
            uniqueEdges.push_back(sortedEdge);
        }
    }

    // Print the unique edges
    cout << "Union of edges in superior triangles and source edge set:\n";
    for (int i = 0; i < uniqueEdges.size(); i++) {
        cout << "(" << uniqueEdges[i].first << ", " << uniqueEdges[i].second << ")\n";
    }
}
//...........................
void graph::optimizedKtrussForUnionEdges() {
    // Clear previous results
    results.clear();

    // Iterate over the edges in the union edges
    for (int i = 0; i < uniqueEdges.size(); ++i) {
        int u = uniqueEdges[i].first;
        int v = uniqueEdges[i].second;

        // Calculate common neighbors of u and v
        int commonCount = 0;
        int idxU = 0, idxV = 0;
        while (idxU < neighborList[u].size() && idxV < neighborList[v].size()) {
            if (neighborList[u][idxU] == neighborList[v][idxV]) {
                ++commonCount;
                ++idxU;
                ++idxV;
            }
            else if (neighborList[u][idxU] < neighborList[v][idxV]) {
                ++idxU;
            }
            else {
                ++idxV;
            }
        }

        // Store the result (u, v, K-truss)
        //vector < int > resultEntry = {u, v, commonCount};
        vector<int> resultEntry;
        resultEntry.push_back(u);
        resultEntry.push_back(v);
        resultEntry.push_back(commonCount);
        results.push_back(resultEntry);
    }

    // Print the results
    cout << "K-truss values for edges in the union:\n";
    for (const auto& res : results) {
        cout << "(" << res[0] << ", " << res[1] << ", K=" << res[2] << ")\n";
    }
}

//...........................
void graph::findKinUnionEdges() {
    // Vector to store results (a, b, K)
    results.clear();

    // Add edges from superior triangles
    for (int i = 0; i < superiorTriangles.size(); i++) {
        vector <int> triangle = superiorTriangles[i];
        vector< pair <int, int> > edges;
        edges.push_back(pair<int, int>(min(triangle[0], triangle[1]), max(triangle[0], triangle[1])));
        edges.push_back(pair<int, int>(min(triangle[1], triangle[2]), max(triangle[1], triangle[2])));
        edges.push_back(pair<int, int>(min(triangle[0], triangle[2]), max(triangle[0], triangle[2])));

        /*vector <pair <int, int>> edges = {
            {min(triangle[0], triangle[1]), max(triangle[0], triangle[1])},
            {min(triangle[1], triangle[2]), max(triangle[1], triangle[2])},
            {min(triangle[0], triangle[2]), max(triangle[0], triangle[2])}
        };*/

        for (auto& edge : edges) {
            if (!edgeExists(uniqueEdges, edge)) {
                uniqueEdges.push_back(edge);
                // Find K for this edge
                int a = edge.first, b = edge.second;
                int commonCount = 0;

                // Count common neighbors of a and b
                for (int x : neighborList[a]) {
                    for (int y : neighborList[b]) {
                        if (x == y) {
                            ++commonCount;
                        }
                    }
                }

                // Store result (a, b, K)
                vector < int> temp1;
                temp1.push_back(a);
                temp1.push_back(b);
                temp1.push_back(commonCount);
                results.push_back(temp1);
                results.push_back(temp1);
            }
        }
    }

    // Add edges from the source edge set
    for (int i = 0; i < sourceEdgeSet.size(); i++) {
        pair<int, int> edge = sourceEdgeSet[i];
        //pair<int, int> sortedEdge = {min(edge.first, edge.second), max(edge.first, edge.second)};
        int smaller = min(edge.first, edge.second);
        int larger = max(edge.first, edge.second);
        pair<int, int> sortedEdge(smaller, larger);

        if (!edgeExists(uniqueEdges, sortedEdge)) {
            uniqueEdges.push_back(sortedEdge);
            // Find K for this edge
            int a = sortedEdge.first, b = sortedEdge.second;
            int commonCount = 0;

            // Count common neighbors of a and b
            for (int x : neighborList[a]) {
                for (int y : neighborList[b]) {
                    if (x == y) {
                        ++commonCount;
                    }
                }
            }

            // Store result (a, b, K)
            vector < int> temp;
            temp.push_back(a);
            temp.push_back(b);
            temp.push_back(commonCount);
            results.push_back(temp);
        }
    }

    // Print the unique edges with their K values
    cout << "Union of edges with common neighbors (a, b, K):\n";
    for (int i = 0; i < results.size(); i++) {
        cout << "(" << results[i][0] << ", " << results[i][1] << ", " << results[i][2] << ")\n";
    }
}

//...........................
void graph::findSourceEdgeSet(int a, int b) {
    vector <pair <int, int> > sourceEdgeSet;

    // Get neighbors of a and b
    vector<int> neighborsA = getNeighbors(a);
    vector<int> neighborsB = getNeighbors(b);

    // Find common neighbors (without sorting)
    vector<int> commonNeighbors;
    for (int neighborA : neighborsA) {
        for (int neighborB : neighborsB) {
            if (neighborA == neighborB) {
                commonNeighbors.push_back(neighborA);
            }
        }
    }

    // Add edges involving common neighbors
    for (int neighbor : commonNeighbors) {
        addEdge(sourceEdgeSet, a, neighbor);
        addEdge(sourceEdgeSet, b, neighbor);
    }

    // Add edge (a, b) itself
    addEdge(sourceEdgeSet, a, b);

    cout << "Source Edge Set for edge (" << a << ", " << b << "):\n";
    for (const auto& edge : sourceEdgeSet) {
        cout << "(" << edge.first << ", " << edge.second << ")\n";
    }
}
//...........................
void graph::setGraph(int totalNode, int totalEdge) {

    this->totalNode = totalNode;
    this->totalEdge = totalEdge;
    csr.resize(totalNode + totalEdge, -1);
    degree.resize(totalNode, -1);
    supports.resize(totalEdge / 2, vector<int>(3, -1));
}
//...........................
void graph::processTriangleEdges(const vector<int>& triangle) {
    // Ensure triangle has 3 nodes
    if (triangle.size() != 3) {
        cout << "Invalid triangle passed to processTriangleEdges.\n";
        return;
    }
    // Extract edges from the triangle
    vector <pair <int, int> > edges;
    edges.push_back(make_pair(triangle[0], triangle[1]));
    edges.push_back(make_pair(triangle[1], triangle[2]));
    edges.push_back(make_pair(triangle[0], triangle[2]));
    // Print or process edges
    cout << "Edges of Triangle (" << triangle[0] << ", " << triangle[1] << ", " << triangle[2] << "): ";
    for (const auto& edge : edges) {
        cout << "(" << edge.first << ", " << edge.second << ") ";
    }
    cout << "\n";

    // Additional processing can go here, if needed
}
//.........................
void graph::addEdgeToNeighborList(int u, int v) {
    auto& neighborsU = neighborList[u];
    if (find(neighborsU.begin(), neighborsU.end(), v) == neighborsU.end()) {
        auto it = neighborsU.begin();
        while (it != neighborsU.end() && *it < v) ++it;
        neighborsU.insert(it, v); // Insert in sorted position
    }

    auto& neighborsV = neighborList[v];
    if (find(neighborsV.begin(), neighborsV.end(), u) == neighborsV.end()) {
        auto it = neighborsV.begin();
        while (it != neighborsV.end() && *it < u) ++it;
        neighborsV.insert(it, u); // Insert in sorted position
    }
}

//..........................
void graph::generateneighborList() {
    neighborList.resize(totalNode);
    for (int i = 0; i < totalNode; ++i) {
        int start = csr[i];
        int end = (i == totalNode - 1) ? csr.size() : csr[i + 1];
        for (int j = start; j < end; ++j) {
            neighborList[i].push_back(csr[j]);
        }
        // Sort neighbors for each node
        //sort(neighborList[i].begin(), neighborList[i].end());
    }
}
//...........................
void graph::findSuperiorTriangles() {
    cout << "Finding Superior Triangles...\n";
    // Vector to store unique triangles
    for (int s = 0; s < supports.size(); ++s) {
        int a = supports[s][0];
        int b = supports[s][1];
        int supportAB = supports[s][2];
        // Truss number for edge (a, b)
        int trussAB = supportAB + 2;
        // Neighbors of 'a' and 'b' from CSR
        int startIndexA = csr[a];
        int endIndexA = (a == totalNode - 1) ? csr.size() : csr[a + 1];
        int startIndexB = csr[b];
        int endIndexB = (b == totalNode - 1) ? csr.size() : csr[b + 1];
        for (int i = startIndexA; i < endIndexA; ++i) {
            int c = csr[i]; // Neighbor of 'a'
            if (c == b) continue;
            // Check if 'c' is also a neighbor of 'b'
            for (int j = startIndexB; j < endIndexB; ++j) {
                if (csr[j] == c) {
                    // Triangle (a, b, c) found
                    // Truss numbers for edges (b, c) and (a, c)
                    int trussBC = 0, trussAC = 0;
                    // Find support for edge (b, c)
                    for (int k = 0; k < supports.size(); ++k) {
                        if ((supports[k][0] == b && supports[k][1] == c) || (supports[k][0] == c && supports[k][1] == b)) {
                            trussBC = supports[k][2] + 2;
                            break;
                        }
                    }
                    // Find support for edge (a, c)
                    for (int k = 0; k < supports.size(); ++k) {
                        if ((supports[k][0] == a && supports[k][1] == c) || (supports[k][0] == c && supports[k][1] == a)) {
                            trussAC = supports[k][2] + 2;
                            break;
                        }
                    }
                    // Check Superior Triangle condition
                    if ((trussAB < trussBC && trussAB < trussAC) ||
                        (trussBC < trussAB && trussBC < trussAC) ||
                        (trussAC < trussAB && trussAC < trussBC)) {
                        // Create a sorted triangle
                        vector<int> triangle;
                        triangle.push_back(a);
                        triangle.push_back(b);
                        triangle.push_back(c);
                        for (int m = 0; m < 2; ++m) {
                            for (int n = m + 1; n < 3; ++n) {
                                if (triangle[m] > triangle[n]) {
                                    int temp = triangle[m];
                                    triangle[m] = triangle[n];
                                    triangle[n] = temp;
                                }
                            }
                        }
                        // Check for uniqueness using simple loop
                        bool isUnique = true;
                        for (int k = 0; k < superiorTriangles.size(); ++k) {
                            if (superiorTriangles[k] == triangle) {
                                isUnique = false;
                                break;
                            }
                        }
                        // Add triangle if unique
                        if (isUnique) {
                            superiorTriangles.push_back(triangle);
                            processTriangleEdges(triangle);
                        }
                    }
                }
            }
        }
    }
    // Print the unique superior triangles
    cout << "Unique Superior Triangles:\n";
    for (int i = 0; i < superiorTriangles.size(); ++i) {
        cout << "(" << superiorTriangles[i][0] << ", " << superiorTriangles[i][1] << ", " << superiorTriangles[i][2] << ")\n";
    }
}
//...........................
void graph::Ktruss() {
    int from, to, a, b, bi, startIndexA, endIndexA, startIndexB, endIndexB, count, s = 0;
    for (int a = 0; a < totalNode; a++) {
        from = csr[a];
        to = csr[a + 1];
        if (a == totalNode - 1)
            to = csr.size();
        for (int bi = from; bi < to; bi++) {//neighbors of a
            b = csr[bi];
            if (a < b)
                continue;
            count = 0;
            //___________________
            startIndexA = csr[a];
            endIndexA = csr[a + 1];
            if (a == totalNode - 1)
                endIndexA = csr.size();
            startIndexB = csr[b];
            endIndexB = csr[b + 1];
            if (b == totalNode - 1) {
                endIndexB = csr.size();
            }
            //___________________
            int i = startIndexA;
            int j = startIndexB;
            do
            {
                if (csr[i] == csr[j]) {
                    count++;
                    i++;
                    j++;
                }
                else if (csr[i] < csr[j]) {
                    i++;
                }
                else {
                    j++;
                }
            } while (i < endIndexA && j < endIndexB);
            //cout<<" A: "<<a<<"B: "<<b<<" Support: " <<count<<endl;

            supports[s][0] = a;
            supports[s][1] = b;
            supports[s][2] = count;
            s++;
        }
    }
}
//...........................
void graph::findCommonNeighbors() {
    for (int a = 0; a < totalNode; ++a) {
        for (int i = 0; i < neighborList[a].size(); ++i) {
            int b = neighborList[a][i];
            if (a < b) {
                // Find common neighbors manually
                vector<int> commonNeighbors;
                const auto& neighborsA = neighborList[a];
                const auto& neighborsB = neighborList[b];
                int idxA = 0, idxB = 0;
                while (idxA < neighborsA.size() && idxB < neighborsB.size()) {
                    if (neighborsA[idxA] == neighborsB[idxB]) {
                        commonNeighbors.push_back(neighborsA[idxA]);
                        ++idxA;
                        ++idxB;
                    }
                    else if (neighborsA[idxA] < neighborsB[idxB]) {
                        ++idxA;
                    }
                    else {
                        ++idxB;
                    }
                }

                int K = commonNeighbors.size();
                //vector <int> temp = { a, b, K };
                vector<int> temp;
                temp.push_back(a);
                temp.push_back(b);
                temp.push_back(K);
                results.push_back(temp);
            }
        }
    }
    for (const auto& row : results) {
        cout << "(" << row[0] << ", " << row[1] << ", " << row[2] << ")\n";
    }
}
//............................
void graph::printSupport() {
    for (int i = 0; i < supports.size(); i++) {
        cout << "node= " << supports[i][0] << "    " << "node= " << supports[i][1] << " support: " << supports[i][2] << endl;
    }
}
//...........................
void printGraph(graph g) {
    cout << endl;
    int from, to;
    for (int i = 0; i < g.totalNode; i++) {
        from = g.csr[i];
        to = g.csr[i + 1];
        if (i == g.totalNode - 1) {
            to = g.csr.size();
        }
        cout << i << ": ";
        for (int j = from; j < to; j++) {
            cout << g.csr[j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    //printDegrees(g.degree);
}
//..........................
void graph::printGraphVector() {
    cout << "Sorted Neighbor List for Each Node:" << endl;
    for (int i = 0; i < neighborList.size(); ++i) {
        cout << "Node " << i << ": ";
        for (int neighbor : neighborList[i]) {
            cout << neighbor << " ";
        }
        cout << endl;
    }
}
//..........................
void readFile(graph& g, string fileName) {
    fstream p;
    int src, dest, totalEdge, totalNode, rowIndex = 0, lastSrc = -1;
    p.open(fileName);
    p >> totalNode >> totalEdge;
    g.setGraph(totalNode, totalEdge);
    int colIndex = totalNode;
    while (!p.eof()) {
        p >> src >> dest;
        if (dest == -1) {
            rowIndex++;
            lastSrc = src;
            g.degree[src] = -1;
            continue;
        }
        else if (src != lastSrc) {
            g.csr[rowIndex] = colIndex;
            rowIndex++;
            lastSrc = src;
        }
        g.csr[colIndex] = dest;
        colIndex++;
        g.degree[src]++;
    }
    p.close();
    // for(int i=0;i<g.csr.size();i++){
     //    cout<<i<< " : " <<g.csr[i]<<endl;
    // }
}
//..........................................
int main() {
    graph g;
    string fileName = "net2.txt";
    readFile(g, fileName);
    printGraph(g);
    g.Ktruss();
    g.printSupport();
    g.generateneighborList();
    g.findSuperiorTriangles();
    //..........................
    int a = 2, b = 8; // Edge to check
    g.addEdgeToNeighborList(a, b);
    g.findSourceEdgeSet(a, b);
    g.findUnionEdges();
    g.findKinUnionEdges();
    g.printGraphVector();
    g.findCommonNeighbors();
    g.optimizedKtrussForUnionEdges();
}