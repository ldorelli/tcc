#include <kuramoto.h>
#include <cstdlib> 
#include <cmath>
#include <util.h>
#include <vector>

KuramotoOscillator::KuramotoOscillator(double w, double k, double phase, int N) :
	 K(k), N(N) {
	this->w = w;	
	this->phase = phase;
	delta = phase;
}

void KuramotoOscillator::applyStep () {
	phase += delta;
	while ( fabs(phase) > 1.0) phase /= 10.0;
	if (phase < 0) phase = -phase;
	delta = w;
}

void KuramotoOscillator::computeStep(KuramotoOscillator& o) {
	delta += (double) K/N * sin(phase - o.phase);
}

void KuramotoOscillator::simpleKuramotoNetwork(
	Graph<KuramotoOscillator>& g, int n, double k) {
	g.adj = std::vector< std::vector<int> > (n);
	for (int i = 0; i < n; ++i) {
		g.nodes.push_back(KuramotoOscillator(normal(0.0, 1.0), k, 
			-1.0 + 2.0*(double) rand()/RAND_MAX, n) ); 
	}
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j)
			if (j != i) g.adj[i].push_back(j);	
	}
}

