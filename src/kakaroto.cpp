#include "igraph.h"
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <string>
#include "util.hpp"

using namespace std;

class Kakaroto{
public:
	igraph_t graph;
	vector< vector< double > > theta;
	vector< double > t;
	vector<double> omega;
	double sigma;

	//nao le o grafo

	Kakaroto () {}

	Kakaroto (vector< double > _theta0, double _t0, vector< double > _omega, double _sigma) {
		theta.resize(_theta0.size());
		for (int i = 0; i < _theta0.size(); i++)
			theta[i].push_back(_theta0[i]);
		t.push_back(_t0);
		omega = _omega;
		sigma = _sigma;
	}

	double f (int curr, vector<double> k, double coef) {
		double sum = 0;
		int i, size, next;
		igraph_vector_t nid;

		igraph_vector_init (&nid, 0);
		igraph_neighbors(&graph, &nid, curr, IGRAPH_IN);

		size = igraph_vector_size(&nid);
		for (i = 0; i < size; i++) {
			if (i == curr)	continue;
			next = (int)VECTOR(nid)[i];
			sum += sin ((theta[next][theta[next].size()-1]+k[next]*coef)-(theta[curr][theta[curr].size()-1]+k[curr]*coef));
		}
		return omega[i]+sigma*sum;
	}
	
	double f2 (int curr, vector<double> k, double coef) {
		double sum = 0;
		int i, size, next;
		return theta[curr][theta[curr].size()-1]+k[curr]*coef;	
	}

	void calc (int n) {
		int i, j, size;
		vector<double> k0, k1, k2, k3, k4;
		double h = 1e-4, ans, aa;
		size = igraph_vcount (&graph);
		k0 = vector<double> (size, 0);
		for (i = 1; i <= n; i++) {
			k1.clear();
			k2.clear();
			k3.clear();
			k4.clear();
			for (j = 0; j < size; j++)
				k1.push_back(f(j, k0, 0));
			for (j = 0; j < size; j++)
				k2.push_back(f(j, k1, h/2));
			for (j = 0; j < size; j++)
				k3.push_back(f(j, k2, h/2));
			for (j = 0; j < size; j++)
				k4.push_back(f(j, k3, h));
			for (j = 0; j < size; j++) {
				ans = theta[j][i-1] + h/6*(k1[j]+2*k2[j]+2*k3[j]+k4[j]);
				aa = fabs(ans)/(2*M_PI);
				if (ans < 0)	ans += aa+2*M_PI;
				else	ans -= aa;
				theta[j].push_back (ans);
			}
		}
	}
	
	void print (void) {
		int i, j, size;
		string s = "th0";
		size = igraph_vcount(&graph);
		for (i = 0; i < size; i++) {
			Util::printRvector(theta[i], s);
			s[2]++;
		}
	}
};

#include "util.hpp"


int main (void) {
	Kakaroto goku;
	vector<double> theta, omega;
	double t0, sigma;
	string fn = "../networks/teste.el";
	srand(time(NULL));
	theta.push_back(0);
	theta.push_back(1);
	omega.push_back(-0.5);
	omega.push_back(0.3);
	t0 = 0;
	sigma = 0.5;
	goku = Kakaroto (theta, t0, omega, sigma);
	cerr << "plol" << endl;
	Util::readGraph (&(goku.graph), 0, 0, fn.c_str());
	cerr << "lala" << endl;
	goku.calc(10000);
		cerr << "lalal" << endl;
	goku.print();
	return 0;
}
