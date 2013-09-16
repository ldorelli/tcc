#include "igraph.h"
#include "stat.hpp"

class Graph {
public:
	static void barabasi_game(igraph_t *graph, int population, int param, bool isDirected) {
		vector<double> p (param, isDirected);
		vector<int> id;
		int i, j, from, to;
		igraph_empty(graph, population, isDirected);
		for (i = 0; i < param; i++) {
			for (j = i+1; j < param; j++) {
				from = j;
				to = i;
				igraph_add_edge (graph, from, to);
				p[to]++;
				if (!isDirected) p[from]++;
			}
		}
		for (; i < population; i++) {
			Statistic::sample(p, param, id);
			p.push_back(isDirected);
			for (j = 0; j < param; j++) {
				from = i;
				to = id[j];
				igraph_add_edge (graph, from, to);
				p[to]++;
				if (!isDirected)	p[from]++;
			}
		}
	}
};
