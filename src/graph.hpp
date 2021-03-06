#include "igraph.h"
#include "stat.hpp"

class Graph {
public:
	static void SF2ER (igraph_t *graph, int N, int param, double alpha, bool isDirected) {

		int i, j, k, m;
		igraph_bool_t areconn;
		vector<int> p (N, 0); //quantidade de arestas PA
		int etot = 0, acc;
		igraph_empty(graph, N, isDirected);
		for (i = 0; i < 2*param+1; i++)
			for (j = i+1; j < 2*param+1; j++) {
				igraph_add_edge (graph, i, j);
				p[j]++;
				if (!isDirected)	p[i]++;
				etot += 1+(!isDirected);
			}
		for (; i < N; i++) {			
			for (m = param; m > 0; m--) {
				if (Statistic::unif() < alpha) { //ER
					do {
						j = rand()%(N-1);
						if (j >= i)	j++;
						igraph_are_connected(graph, i, j, &areconn);
						//printf ("er %d\n", (int)areconn);
					} while (areconn == 1);
					igraph_add_edge (graph, i, j);
				} else { //SF
					do {
						k = rand()%(etot-p[i]);
						//printf ("%d %d %d\n", k, etot, p[i]);
						acc = 0;
						for (j = 0; j < N; j++) {
							if (i == j)	continue;
							if (acc+p[j] >= k)	break;
							acc+=p[j];
						}
						igraph_are_connected(graph, i, j, &areconn);
						//printf ("sf %d %d %d\n",k, (int)areconn, etot);
					} while (areconn == 1);
					igraph_add_edge (graph, i, j);
					p[j]++;
					if (!isDirected)	p[i]++;
					etot += 1+(!isDirected);
				}
			}
		}
	}
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
	}
};
