
#include "Graph.hpp"
#include <vector>
#include <iostream>
using namespace std;//there was no error but i added
namespace ariel{
    Graph::Graph(): mat(), isDirected(){}
    
    bool directed;
    //isSymetric = is directed 
    bool isSymetric(const vector<vector<int>>& ee) {
        using vector_size_type = typename vector<vector<int>>::size_type;
        for (vector_size_type i = 0; i < ee.size(); ++i) {
            for (vector_size_type j = 0; j < ee.size(); ++j) {
                if (ee[i][j] != ee[j][i]) return false;
            }
        }
        return true;
    }

    void Graph::loadGraph(const vector<vector<int>>& gra3){
        if (gra3.size()!=gra3[0].size()){
            cerr << "ERROR loading graph. Entered matrix isn't square. please try again. \n";//Checking N*N
        }
        else {
            bool symetric = (isSymetric(gra3)==true) ? true : false;
            if (symetric ==false) directed = false;
            else {//matrix is symetric
                char choice;
                cout << "I've noticed the matrix is symmertic. Is your graph directed? y\n" << endl;
                cin >> choice;
                directed = (choice == 'y') ? true : false;
            }
            isDirected = directed;
            mat = gra3;//actual load of matrix
        }
    }

    void Graph::printGraph(){
        int n = Graph::mat.size();
        int num_edge = countEdges(mat, isDirected);
        cout << "Graph with " << n << "vertices and " << num_edge << " edges" << endl;//check with gpt
    }

    int countEdges(vector<vector<int>> graph, bool directed) {
        int numEdges = 0;
        int numRows = graph.size();
        int numCols = graph[0].size();

        for (unsigned int i = 0; i < numRows; ++i) {
            for (unsigned int j = 0; j < numCols; ++j) {
                if (graph[i][j] != 0) {
                    if (!directed && j <= i) continue;// For undirected graphs, count only one direction
                    numEdges++;
                }
            }
        }
        // If undirected, divide by 2 to avoid double counting
        //if (!directed) {
        //    numEdges /= 2;
        //}
        return numEdges;
    }


}