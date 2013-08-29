#include "igraph.h"
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <SFML/Graphics.hpp>
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
		double h = 1e-1, ans, aa;
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
				ans = theta[j][i-1] + (h/6.0)*(k1[j]+2*k2[j]+2*k3[j]+k4[j]);
				aa = floor( fabs(ans)/(2*M_PI) );
				if (ans < 0)	ans += aa+2*M_PI;
				else	ans -= aa*2*M_PI;
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
	theta.push_back(2.33);
	omega.push_back(-0.5);
	omega.push_back(1.20);
	omega.push_back(0.3);
	t0 = 0;
	sigma = 0.01;
	goku = Kakaroto (theta, t0, omega, sigma);
	Util::readGraph (&(goku.graph), 0, 0, fn.c_str());
	goku.calc(100000);
	// goku.print();

	sf::RenderWindow window(sf::VideoMode(800, 600), "My window");
	sf::View view(
		sf::Vector2f(0.0, 0.0), 
		sf::Vector2f(200,150) );
	
	window.setView(view);

	// run the program as long as the window is open    
	for (int i = 0; i < 100000; ++i) {
    	sf::Event event;
    	while (window.pollEvent(event))
    	{
        	// "close requested" event: we close the window
        	if (event.type == sf::Event::Closed)
            	window.close();
    	}
    	// clear the window with black color
    	window.clear(sf::Color::Black);
    	int size = igraph_vcount(&goku.graph);
    	for (int j = 0; j < size; ++j) {
    		double tt = goku.theta[j][i];
    		sf::CircleShape sp(2);
			sp.setPosition(-20+5*j, 0); 
			sp.setFillColor( sf::Color(255*tt/(2*M_PI),
				255*tt/(2*M_PI), 255*tt/(2*M_PI)) );
			window.draw(sp);
    	}
    	window.display();
    	sf::sleep( sf::seconds(0.00001) );	
    }


	return 0;
}