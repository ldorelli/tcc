#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp> 
#include <vector>
#include <iostream>

// Oscillator indica o tipo do nó
template <class Oscillator, class Phase>
class Graph {
public:
	std::vector<Oscillator>	nodes;
	std::vector< std::vector<int> > adj;
	// Applies all changes 
	void update();

private:
};

template <class Oscillator, class Phase>
void Graph<Oscillator, Phase>::update() {
	Phase m, M;
	for (int i = 0; i < nodes.size(); ++i) {
		for (int j = 0; j < adj[i].size(); ++j)
			nodes[i].computeStep(nodes[adj[i][j]]);
	}
	m = M = nodes[0].phase;
	for (int i = 0; i < nodes.size(); ++i) {
		nodes[i].applyStep();
		std::cout << nodes[i].phase << " ";
		m = min (m, nodes[i].phase);
		M = max (M, nodes[i].phase);
	}
	std::cout << std::endl;
//	for (int i = 0; i < nodes.size(); ++i)
//		nodes[i].normalize(m,M);
}

#endif 
