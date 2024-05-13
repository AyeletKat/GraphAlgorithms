#ifndef Algorithms_HPP
#define Algorithms_HPP
#include <string>
#include "Graph.hpp"


namespace ariel {
    //class Graph;
    class Algorithms{
        public:
            static int isConnected(const ariel::Graph& g1);
            static int shortestPath(const ariel::Graph& g2, const int start, const int end);//returns the SP or vertices or -1.
            static int isContainsCycle(const ariel::Graph& g3);//prings cycle or 0
            static int isBipartite(const ariel::Graph& g4);//returns two groups of vertices or 0
            static int negativeCycle(const ariel::Graph& g5);//find negative cycle or print if there is not
        private:
            Algorithms(){}
    };
};
#endif