#include "igraph.h"
#include <cmath>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <queue>
#include <set>
#include "util.hpp"

using namespace std;

int main (int argc, char* argv[]) {
	char *p;
	igraph_t graph;
	if (argc < 3) {
		fprintf (stderr, "Usage: %s <1-filename> <2-properties: 0.<k>; 1.C; 2.<k^2>/<k>; 3.l; 4.connectivity> <3-degree distribution filename>\n", argv[0]);
		return 1;
	}
	Util::readGraph(&graph, argv[1]);
	p = argv[2];
	if (p[0] == '1') {
		int e = (int)igraph_ecount(&graph);
		int n = (int)igraph_vcount(&graph);
		double k = 2.0*e/n;
		cout << "<k> = "<< k << endl;
	} if (p[1] == '1') {
		igraph_real_t clust;
		igraph_transitivity_undirected(&graph, &clust, IGRAPH_TRANSITIVITY_ZERO);
		cout << "C = " << (double) clust << endl;
	} if (p[2] == '1') {
		igraph_vector_t v_deg;
		vector<double> vDeg;
		map<double, double> dDeg;
		igraph_vector_init(&v_deg, 0);
		igraph_degree(&graph, &v_deg, igraph_vss_all(), IGRAPH_ALL, 0);
		Util::makeVector(&v_deg, vDeg);
		Statistic::getDistribution(vDeg, dDeg);
		double k = Statistic::getMoment(1, dDeg);
		double k2 = Statistic::getMoment(2, dDeg);
		double h = k2/k;
		cout << "h = " << h << endl;
	} if (p[3] == '1') {
		igraph_real_t avgsp;
		igraph_average_path_length(&graph, &avgsp, IGRAPH_UNDIRECTED, 0);
		cout << "l = " << avgsp << endl;
	} if (p[4] == '1') {
		igraph_bool_t conn;
		igraph_is_connected(&graph, &conn, IGRAPH_STRONG);
		cout << "isCon = " << (int)conn << endl;
	}
	if (argc > 3) {
		igraph_vector_t v_deg;
		vector<double> vDeg;
		igraph_vector_init(&v_deg, 0);
		igraph_degree(&graph, &v_deg, igraph_vss_all(), IGRAPH_ALL, 0);
		Util::makeVector(&v_deg, vDeg);
		Util::printRvector(string(argv[4]), vDeg, "degree");
	}
	return 0;
}
