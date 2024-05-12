# ifndef Graph_HPP
# define Graph_HPP

# include <vector>

using namespace std;

namespace ariel {

    class Graph{
        //(private)

    public: 
        bool isDirected;
        vector<vector<int>> mat;
        Graph();
        Graph(const vector<vector<int>>& g1){//inline
            mat = g1;
        }
        void loadGraph(const vector<vector<int>>& gra3);//{//outline//mat = gra3;}
        void printGraph();

        int countEdges(const vector<vector<int>>& graph, bool directed);
        void setDirected(bool b){
            isDirected = b;
        }
    };

};
# endif