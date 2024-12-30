#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>//Random 
#include <ctime>// reset kardan 
using namespace std;

class graph {
public:
    int nodes;
    int edges;
   
    vector< vector <bool> > adjMatrix; // Adjacency matrix for reading graph and put there
    vector< vector <short> > Ktruss; // store k-truss values 
    vector< vector <bool> > sourceEdgesSets;
    vector< vector <bool> > superiorTriangleEdges;

    void ReadFile(string);
    void ComputeKTruss();
    void PrintKTrussMatrix();
    void FindSuperiorTriangles();
    void UpdateGraph();
    void PrintAdjMatrix();
    void unionSet();
    void PrintSuperiorSource();
    void UpdateKTruss(short i,short j);
    void FindSourceSet(short,short);
};
//......................................
void graph::ReadFile(string fileName) {
    ifstream inputFile(fileName);
    if (!inputFile){
        cout << "Error opening file!" << endl;
        return;
    }
    inputFile >> nodes >> edges;
    adjMatrix.resize(nodes, vector<bool>(nodes, false));
    Ktruss.resize(nodes, vector<short>(nodes, 0));
    sourceEdgesSets.resize(nodes, vector<bool>(nodes, false));
    superiorTriangleEdges.resize(nodes, vector<bool>(nodes, false));
    int u, v;
    while (inputFile >> u >> v) {
        adjMatrix[u][v] = true;
        adjMatrix[v][u] = true;
    }
    inputFile.close();
}
//.................................
void graph::UpdateGraph(){
    int i,j;
    i=rand()%nodes;
    j=rand()%nodes;
    //i=4;
    //j=5;
    adjMatrix[i][j]=!adjMatrix[i][j];//
    adjMatrix[j][i]=!adjMatrix[j][i];//bayad not konim ta edge tekrari nashe
    cout<<"\n Updating Edge: ("<<i<<" , "<<j<<" )"<<endl;
    if(adjMatrix[i][j]){
        UpdateKTruss(i,j);//just calculate the k truss i j 
    }
    else{
        Ktruss[i][j]=0;
        Ktruss[j][i]=0;
    }
    FindSourceSet(i,j);
}
//................................
void graph::PrintAdjMatrix() {
    cout << "Adjacency Matrix:" << endl;
    for (int i = 0; i < nodes; ++i) {
        for (int j = 0; j < nodes; ++j) {
            cout << adjMatrix[i][j] << " ";
        }
        cout << endl;
    }
}
//.................................
void graph::UpdateKTruss(short i,short j) {
    int count;
    if (adjMatrix[i][j]) {
        count = 0;
        for (int k = 0; k < nodes; ++k) {
            if (adjMatrix[i][k] && adjMatrix[j][k]) {
                count++;
                }
            }
            Ktruss[i][j] = count + 2;
            Ktruss[j][i] = count + 2;
            }
}
//.................................
void graph::ComputeKTruss() {
    int count;
    for (int i = 0; i < nodes; ++i) {
        for (int j = i + 1; j < nodes; ++j) {//just check upper place of main diagonal to avoid duplicate 
            if (adjMatrix[i][j]) {//if there is an edge
                count = 0;
                for (int k = 0; k < nodes; ++k) {
                    if (adjMatrix[i][k] && adjMatrix[j][k]) {//find common neighbors
                        count++;
                    }
                }
                Ktruss[i][j] = count + 2;
                Ktruss[j][i] = count + 2;
            }
        }
    }
}
//........................................
void graph::PrintKTrussMatrix() {
    cout << "K-Truss Adjacency Matrix:" << endl;
    for (int i = 0; i < nodes; ++i) {
        for (int j = 0; j < nodes; ++j) {
            cout << Ktruss[i][j] << " ";
        }
        cout << endl;
    }
}
//..........................................
void graph::FindSuperiorTriangles() {
    cout << "Superior Triangles:" << endl;
    for (int i = 0; i < nodes; ++i) {
        for (int j = i + 1; j < nodes; ++j) {
            if (Ktruss[i][j]>2) {//if i , j member of at least 1 triangle 
                for (int k = j + 1; k < nodes; ++k) {// balaye ghotr asli
                    if (adjMatrix[i][k] && adjMatrix[j][k]) {// az i be k edge vojood dare va j be k , common neighbor vojood dare 
                        // Check for superior triangle
                        short edge1 = Ktruss[i][j];// k truss i j 
                        short edge2 = Ktruss[i][k];
                        short edge3 = Ktruss[j][k];
                        if ((edge1 > edge3 && edge2 > edge3 && edge1==edge2) ||
                            (edge1 > edge2 && edge3 > edge2 && edge1==edge3) ||
                            (edge2 > edge1 && edge3 > edge1 && edge2==edge3)) {
                            superiorTriangleEdges[i][j]=true; //mark the edges as superior triangle and print
                            superiorTriangleEdges[j][i]=true;
                            cout << "Triangle: (" << i << ", " << j  <<")" << endl;
                            superiorTriangleEdges[i][k]=true;
                            superiorTriangleEdges[k][i]=true;
                            cout << "Triangle: (" << i << ", " << k << ")" << endl;
                            superiorTriangleEdges[j][k]=true;
                            superiorTriangleEdges[k][j]=true;
                            cout << "Triangle: (" <<  j << ", " << k << ")" << endl;
                        }
                    }
                }
            }
        }
    }
}
//............................
void graph::PrintSuperiorSource() {
    cout << "Source Edges Set:" << endl;
    for (int i = 0; i < nodes; ++i) {
        for (int j = 0; j < nodes; ++j) {
            cout << sourceEdgesSets[i][j] << " ";
        }
        cout << endl;
    }
    cout << "superior Triangle Edges:" << endl;
    for (int i = 0; i < nodes; ++i) {
        for (int j = 0; j < nodes; ++j) {
            cout << superiorTriangleEdges[i][j] << " ";
        }
        cout << endl;
    }
}
//............................
void graph::unionSet(){
    for(int i=0;i<nodes;i++){
        for(int j=i+1;j<nodes;j++){
            if(sourceEdgesSets[i][j]|| superiorTriangleEdges[i][j]){//if each one would be true ,it prints its k truss
                cout<<"( "<<i<<", "<<j<<" ) = "<<Ktruss[i][j]<<endl;
            }
        }
    }
}
//............................
void graph::FindSourceSet(short a,short b){
    for (int i = 0; i < sourceEdgesSets.size(); ++i) {//Before starting,it clears all previous entries in the sourceEdgesSets matrix by setting all to false.
    fill(sourceEdgesSets[i].begin(), sourceEdgesSets[i].end(), false);
    }
    int lastK;
   //cout<<"SourceEdge Set"<<endl;
    for(int i=0;i<nodes;i++){
        if(adjMatrix[a][i] && adjMatrix[b][i]){
            //updating edge(a,i) 
            lastK=Ktruss[a][i];
            UpdateKTruss(a,i);
            if(lastK!=Ktruss[a][i]){//taht tasir gharar gerefte
                sourceEdgesSets[a][i]=true;
                sourceEdgesSets[i][a]=true;
                //cout << " (" <<  a << ", " << i << ")" << endl;
            }
            //updating edge(b,i) 
            lastK=Ktruss[b][i];
            UpdateKTruss(b,i);
            if(lastK!=Ktruss[b][i]){//k truss taht tasir gharar gerefte
                sourceEdgesSets[i][b]=true;
                sourceEdgesSets[b][i]=true;
            //cout << " (" <<  b << ", " << i << ")" << endl;
            }
        }//end of for
    sourceEdgesSets[a][b]=true;
    }
}
//............................
int main() {
    srand(time(0));//generating random edges
    graph g;
    g.ReadFile("network6.txt");
    g.ComputeKTruss();
    //g.PrintKTrussMatrix();
    g.FindSuperiorTriangles();
    g.UpdateGraph();
    //g.PrintKTrussMatrix();
    //g.PrintSuperiorSource();
    g.unionSet();
    return 0;
}