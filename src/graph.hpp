#include "igraph.h"
#include "stat.hpp"

class Graph {
public:
	static void barabasi_game(igraph_t *graph, int population, int param, bool isDirected) {
		vector<double> p (2*param+1, isDirected);
		vector<int> id;
		int i, j, from, to;
		igraph_empty(graph, population, isDirected);
		for (i = 0; i < 2*param+1; i++) {
			for (j = i+1; j < 2*param+1; j++) {
				from = j;
				to = i;
				igraph_add_edge (graph, from, to);
				p[to]++;
				if (!isDirected) p[from]++;
			}
		}
		int k = 0;
		for (i = 0; i < 2*param+1; i++)	k += p[i];
		fprintf (stderr, " awawa %lf\n", k/(double)(2*param+1));
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
		k = 0;
		for (i = 0; i < population; i++)	k += p[i];
		fprintf (stderr, " awawa %lf\n", k/(double)population);
	}
};
