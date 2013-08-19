#ifndef __KURAMOTO_H__
#define __KURAMOTO_H__

#include <graph.h>

class KuramotoOscillator {
public:
	KuramotoOscillator(double w, double K, double phase, int N);
	void computeStep(KuramotoOscillator &o);
	void applyStep();
	static void simpleKuramotoNetwork(
		Graph<KuramotoOscillator>& g, int n, double k);
	
	double K, w;
	double phase; 
	double delta;
	int N;
};

#endif 	
