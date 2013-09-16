#include "igraph.h"
#include <cmath>
#include <vector>
#include <fstream>
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
	vector<double> omega, R;
	double sigma, step;
	int np;


	Kakaroto () {}

	Kakaroto (vector< double > _theta0, double _t0, vector< double > _omega, double _sigma, int _np, double _step) {
		theta.resize(_theta0.size());
		for (int i = 0; i < _theta0.size(); i++)
			theta[i].push_back(_theta0[i]);
		t.push_back(_t0);
		omega = _omega;
		sigma = _sigma;
		np = _np;
		step = _step;
	}

	Kakaroto (string fn) {
		string gr = fn, conf, line;
		double _theta, _omega, _t0, _sigma, _step;
		int size, _np, plo;
		ifstream file;
		gr.append(".gr");
		file.open(gr.c_str());
		getline(file, line);
		Util::readGraph (&graph, line.c_str());
		getline(file, conf);
		cout << conf << endl;
		file.close();
		file.open(conf.c_str());
		file >> _t0 >> _sigma >> _np >> _step;
		t.push_back(_t0);
		sigma = _sigma;
		np = _np;
		step = _step;
		size = igraph_vcount(&graph);
		while (file >> _theta >> _omega && theta.size() < size) {
			theta.push_back(vector<double> ());
			theta[theta.size()-1].push_back(_theta);
			omega.push_back(_omega);
		}
		cerr << size << theta.size() << endl;
		file.close();
		if (theta.size() != size)
			throw -7;
	}

	double f (int curr, vector<double> k, double coef) {
		double sum = 0;
		int i, n, size, next;
		igraph_vector_t nid;

		n = igraph_vcount(&graph);
		igraph_vector_init (&nid, 0);
		igraph_neighbors(&graph, &nid, curr, IGRAPH_IN);

		size = igraph_vector_size(&nid);

		//pacemaker
		if (curr >= n - np)	return omega[curr];

		for (i = 0; i < size; i++) {
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
		double  ans, aa;
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
				k2.push_back(f(j, k1, step/2));
			for (j = 0; j < size; j++)
				k3.push_back(f(j, k2, step/2));
			for (j = 0; j < size; j++)
				k4.push_back(f(j, k3, step));
			for (j = 0; j < size; j++) {
				ans = theta[j][i-1] + (step/6.0)*(k1[j]+2*k2[j]+2*k3[j]+k4[j]);
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

	void calcR() {
		for (int i = 0; i < theta[0].size(); ++i) {
			double r1 = 0.0;
			double r2 = 0.0;
//			r = (1/N^2)*([sum i=1^N cos(theta)]^2 + [sum i=1^N sin(theta)]^2)  23:20
			for (int j = 0; j < theta.size(); ++j) {
				r1 += cos(theta[j][i]), r2 += sin(theta[j][i]);
			}
			double r = r1*r1 + r2*r2;
		//	cout << r << endl;
			R.push_back(sqrt(r)/theta.size());
		}
	}

	void draw_graph (void) {
		sf::RenderWindow window(sf::VideoMode(800, 600), "My window");
		sf::View view(
			sf::Vector2f(0.0, 0.0), 
			sf::Vector2f(200,150) );
		window.setView(view);
		for (int i = 0; i < theta[0].size() && window.isOpen(); ++i) {
			
			sf::Event event;
			double rho = 40;
			
			while (window.pollEvent(event))
			{
				// "close requested" event: we close the window
				if (event.type == sf::Event::Closed)
					window.close();
			}
			window.clear(sf::Color::Black);
			int size = igraph_vcount(&graph);
			
			double angle = 0.0;
			double step = 2*M_PI/(double)size;
			
			vector<double> x(size), y(size);

			for (int j = 0; j < size; ++j) {
				x[j] = rho * cos(angle);
				y[j] = rho * sin(angle);
			}

			for (i = 0; i < size; i++) {
				igraph_vector_t nid;
				igraph_vector_init (&nid, 0);
				igraph_neighbors(&graph, &nid, i, IGRAPH_IN);				
				int adj_sz = igraph_vector_size(&nid);
				for (int j = 0; j < adj_sz; ++j) {

				}
			}

			for (int j = 0; j < size-np; ++j) {
				double tt = theta[j][i];
				sf::CircleShape sp(2);
				sp.setPosition(rho*cos(angle), rho*sin(angle)); 
				sp.setFillColor( sf::Color(255*(j+1)/(double)(size+1),
					255*(j+1)/(double)(size+1), 255*(j+1)/(double)(size+1)) );
				//sp.setFillColor (sf::Color(255, 255, 255));
				window.draw(sp);
			}
		}		
	}

	void draw (string wname) {
		double rho = 10;
		
		sf::RenderWindow window(sf::VideoMode(800, 600), wname.c_str());
		sf::View view(
			sf::Vector2f(0.0, 0.0), 
			sf::Vector2f(200,150) );
		
		window.setView(view);

		// run the program as long as the window is open    
		for (int i = 0; i < theta[0].size() && window.isOpen(); ++i) {
			sf::Event event;
			while (window.pollEvent(event))
			{
				// "close requested" event: we close the window
				if (event.type == sf::Event::Closed)
					window.close();
			}
			// clear the window with black color
			window.clear(sf::Color::Black);
			int size = igraph_vcount(&graph);
			int j;
			for (j = 0; j < size-np; ++j) {
				double tt = theta[j][i];
				sf::CircleShape sp(2);
				sp.setPosition(rho*cos(tt), rho*sin(tt)); 
				sp.setFillColor( sf::Color(255*(j+1)/(double)(size+1),
					255*(j+1)/(double)(size+1), 255*(j+1)/(double)(size+1)) );
				//sp.setFillColor (sf::Color(255, 255, 255));
				window.draw(sp);
			}
			for (; j < size; ++j) {
				double tt = theta[j][i];
				sf::CircleShape sp(2);
				sp.setPosition(rho*cos(tt), rho*sin(tt)); 
				sp.setFillColor( sf::Color(.0,
					.0, 255*(j+1)/(double)(size+1)) );
				//sp.setFillColor (sf::Color(255, 255, 255));
				window.draw(sp);
			}
		//	printf("R: %.5lf\n", R[i]);
			window.display();
			sf::sleep( sf::seconds(0.00001) );	
		}
	}
};

#include "util.hpp"


int main (int argc, char* argv[]) {
	Kakaroto goku;
	vector<double> theta, omega;
	double t0, sigma;
	string fn = "../networks/";
	if (argc > 1)	fn.append(string(argv[1]));
	else	fn.append("plo");
	goku = Kakaroto(fn);
	goku.calc(100000);
	goku.calcR();
	cout << goku.R.back() << endl;
	goku.draw(fn);

	return 0;
}
