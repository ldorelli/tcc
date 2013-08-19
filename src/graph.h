#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp> 
#include <vector>

// Oscillator indica o tipo do nó
template <class Oscillator>
class Graph {
public:
	std::vector<Oscillator>	nodes;
	std::vector< std::vector<int> > adj;
	// Applies all changes 
	void update();

private:
};

template <class Oscillator>
void Graph<Oscillator>::update() {
	for (int i = 0; i < nodes.size(); ++i) {
		for (int j = 0; j < adj[i].size(); ++j)
			nodes[i].computeStep(nodes[adj[i][j]]);
	}
	for (int i = 0; i < nodes.size(); ++i) {
		nodes[i].applyStep();
	}
}

#endif 
