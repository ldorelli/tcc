#ifndef UTIL_H
#define UTIL_H
#include <stdio.h>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <igraph.h>
#include <stdlib.h>
#include "graph.hpp"
using namespace std;

class Util {
public:
	static void readGraph (igraph_t* graph, const char* param) {
		FILE* instream;
		cerr << param << endl;
		if ((instream = fopen (param, "r")) == NULL) {
			fprintf (stderr, "Cannot open file %s\n", param);
			throw 1;
		}
		igraph_read_graph_edgelist(graph, instream, 0, IGRAPH_UNDIRECTED);
		fclose(instream);
	}

	static void genGraph (int POPULATION, int type, char* param, char* fname) {
		int g_m;
		double g_p, alpha, beta;
		char dfn[1024];
		FILE* outstream;
		igraph_t graph;
		switch(type) {
			case 1:
				sscanf (param, "%d", &g_m);
				cerr << "BA " << g_m << endl;
				Graph::barabasi_game(&graph, POPULATION, g_m, false);
				break;
			case 2:
				sscanf (param, "%lf", &g_p);
				cerr << "ER " << g_p << endl;
				igraph_erdos_renyi_game(&graph, IGRAPH_ERDOS_RENYI_GNP, POPULATION, g_p, IGRAPH_UNDIRECTED, false);
				break;
			case 3:
				sscanf (param, "%d,%lf", &g_m, &g_p);
				cerr << "WS " << g_m << " " << g_p << endl;
				igraph_watts_strogatz_game(&graph, 1, POPULATION, g_m, g_p);
				break;
			case 4:
				sscanf (param, "%d,%lf", &g_m, &g_p);
				cerr << "NLBA " << g_m << " " << g_p << endl;
				igraph_nonlinear_barabasi_game(&graph, POPULATION, g_p, g_m, NULL, true, 0.01, IGRAPH_UNDIRECTED);
				break;
			case 5:
				sscanf (param, "%d,%lf", &g_m, &g_p);
				cerr << "SF2ER " << g_m << " " << g_p << endl;
				Graph::SF2ER(&graph, POPULATION, g_m, g_p, IGRAPH_UNDIRECTED);
				break;
			case 6:
				sscanf (param, "%lf,%lf,%s", &alpha, &beta, dfn);
				cerr << "Waxman " << alpha << " " << beta << " " << dfn << endl;
				Graph::Waxman(&graph, POPULATION, alpha, beta, string(dfn), IGRAPH_UNDIRECTED);
				break;
			case 7:
				sscanf (param , "%lf,%lf", &alpha, &beta);
				cerr << "GeoSF(nao implementado) " << alpha << " " << beta << endl;
				break;
			default:
				fprintf (stderr, "Type %d not defined.\n", type);
				return;
		}
		if ((outstream = fopen (fname, "w")) == NULL) {
			fprintf (stderr, "Cannot open file %s\n", fname);
			throw 1;
		}
		int e, n;
		e = (int)igraph_ecount(&graph);
		n = (int)igraph_vcount(&graph);
		fprintf (stderr, "Writing to %s N = %d <k> = %lf\n", fname, n, 2*e/(double)n);
		igraph_write_graph_edgelist(&graph, outstream);
		fclose(outstream);
	}
	
	
	template <typename T>
	static void printRvector (ostream &f, const vector<T> &v, const string & varname) {
		int i = 0;
		f << varname << " = c(";
		for (i = 0; i < v.size()-1; i++)	f << v[i] << ", ";
		f << v[i] << ")" << endl;
	}
	
	template <typename T>
	static void printRvector (string filename, const vector<T> &v, const string &varname) {
		ofstream file;
		file.open(filename.c_str());
		if (!file.good()) {
			cerr << "Error: Could not open " << filename << endl;
			return;
		}
		printRvector(file, v, varname);
		file.close();
	}

	template <typename T>
	static void printRvector (const vector<T> &v, const string & varname) {
		int i;
	
		cout << varname << " = c(";
		for (i = 0; i < v.size()-1; i++)	cout << v[i] << ", ";
		cout << v[i] << ")" << endl;
	}
	template <typename T>
	static void makeVector (igraph_vector_t *iv, vector<T> &v) {
		int i, len = igraph_vector_size(iv);
		v.resize(len);
		for (i = 0; i < len; i++)	v[i] = (T)VECTOR(*iv)[i];
	}	
	
	template <typename F, typename S>
	static void printRmap (map<F, S> &m, const string & varname) {
		typename map<F, S>::iterator it;
	
		cout << varname << " = matrix(c(";
		for (it = m.begin(); it != m.end(); it++)	{
			if (it != m.begin())	cout << ", ";
			cout << it->first << ", " << it->second;
		}
		cout << "), nrow=2)" << endl;
	}
	template <typename F, typename S>
	static void printRmap (ostream &file, map<F, S> &m, const string & varname) {
		typename map<F, S>::iterator it;
	
		file << varname << " = matrix(c(";
		for (it = m.begin(); it != m.end(); it++)	{
			if (it != m.begin())	file << ", ";
			file << it->first << ", " << it->second;
		}
		file << "), nrow=2)" << endl;
	}
	template <typename T>
	static void printRmatrixAsVector (const vector<vector<T> > &M, const string & varname) {
		int i, j;
		
		cout << varname << " = c(";
		for (i = 0; i < M.size(); i++) 
			for (j = 0; j < M[i].size(); j++)
				if (i != M.size()-1 || j != M[i].size()-1) cout << M[i][j] << ", ";
		i--; j--;
		cout << M[j][j] << ")" << endl;
	}
	
	
	static double getAssortativity (igraph_t * graph) {
		int ne, nv, i, j;
		double ki, kj, a, b, c;
		igraph_vector_t v_degree, v_sqdegree;
		igraph_bool_t Aij;

		igraph_vector_init (&v_degree, 0);
		igraph_degree (graph, &v_degree, igraph_vss_all(), IGRAPH_ALL, 0);
		ne = igraph_ecount(graph);
		nv = igraph_vcount(graph);
		a = b = c = 0;
		for (i = 0; i < nv; i++) {
			ki = VECTOR(v_degree)[i];
			for (j = i+1; j < nv; j++) {
				igraph_are_connected(graph, i,j, &Aij);
				if (Aij) {
					kj = VECTOR(v_degree)[j];
					a += ki*kj;
					b += ki + kj;
					c += ki*ki + kj*kj;
				}
			}
		}

		igraph_vector_destroy(&v_degree);
		return (a/ne-0.25*b*b/(ne*ne))/(0.5/ne*b-0.25*c*c/(ne*ne));
	}
	
	static string appendNumber (string s, int n) {
		char x[100];
		sprintf (x, "%d", n);
		s.append(x);
		return s;
	}

	static void permutation (vector<int> &v) {
		int i, curr;
		
		for (i = 0; i < v.size(); i++) {
			curr = rand()%(v.size()-i);
			swap(v[curr], v[v.size()-1-i]);	
		}
	}
};
#endif
