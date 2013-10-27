#ifndef KAKAROTO_H
#define KAKAROTO_H
#include "igraph.h"
#include <cmath>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <queue>
#include <set>
#include <sstream>
#include <SFML/Graphics.hpp>
#include "util.hpp"

using namespace std;

sf::Color RGB_from_freq(double w) {
	if (fabs (w) > M_PI) {
		printf ("THROW %lf\n", w);
		throw -1;
	}
	w += M_PI;
	w = (w/(2*M_PI)) * 400;
	w = 380 + w;
	double R, G, B;
	if (w >= 380 and w < 440) {
    	R = -(w - 440.) / (440. - 380.);
    	G = 0.0;
    	B = 1.0;
	} else if  (w >= 440 and w < 490) {
    	R = 0.0;
    	G = (w - 440.) / (490. - 440.);
    	B = 1.0;
    } else if (w >= 490 and w < 510) {
    	R = 0.0;
    	G = 1.0;
    	B = -(w - 510.) / (510. - 490.)	;
    } else if (w >= 510 and w < 580) {
    	R = (w - 510.) / (580. - 510.);
    	G = 1.0;
    	B = 0.0;
    } else if (w >= 580 and w < 645) {
    	R = 1.0;
    	G = -(w - 645.) / (645. - 580.);
    	B = 0.0;
    } else if (w >= 645 and w <= 780) {
    	R = 1.0;
    	G = 0.0;
    	B = 0.0;
    } else {
    	R = 0.0;
    	G = 0.0;
    	B = 0.0;
    }
    // cout << R << " " << G << " " << B << endl;
    return sf::Color(255*R, 255*G, 255*B);
}

double dif (double a1, double a0) {
	if (a1 > a0)	return min (a1-a0, a0+2*M_PI-a1);
	else	return min (a0-a1, a1+2*M_PI-a1);
}

class Kakaroto{
public:
	igraph_t graph;
	vector< vector< double > > theta, freq;
	vector< double > t;
	vector<double> omega, R, ang, R1, R2; 
	vector<vector<double> > dist2;
	double step;
	int np;
	vector<int> nivel;
	vector<vector<int> > pnivel;
	vector<bool> isPacemaker, isContrarian;
	int delay;
	int stype;
	double phaseDif;

	Kakaroto () {}

	Kakaroto (vector< double > _theta0, double _t0, 
		vector< double > _omega, double _sigma, 
		vector<bool> _isPacemaker, double _step) 
	{
		theta.resize(_theta0.size());
		freq.resize (_theta0.size());
		for (int i = 0; i < _theta0.size(); i++) {
			theta[i].push_back(_theta0[i]);
			freq[i].push_back(_omega[i]);
		}
		t.push_back(_t0);
		omega = _omega;
		dist2.resize(theta.size());
		for (int i = 0; i < theta.size(); i++) {
			dist2[i].resize(theta.size());
			for (int j = 0; j < theta.size(); j++) {
				dist2[i][j] = 1/_sigma;
			}
		}
		isPacemaker = _isPacemaker;
		step = _step;
	}

	Kakaroto (string fn, double _sigma, double _step, int _delay = 0, int stype = 0, int contrarians = 0) {
		string gr = fn, conf, line;
		this->stype = stype;
		double _theta, _omega, _t0;
		int size, _np, plo;
		ifstream file;
		//leitura dos arquivos
		gr.append(".gr");
		file.open(gr.c_str());
		//leitura do grafo
		getline(file, line);
		Util::readGraph (&graph, line.c_str());
		getline(file, conf);
		cerr << conf << endl;
		file.close();
		//leitura da configuracao
		file.open(conf.c_str());
		file >> _t0 >>_np;
		t.push_back(_t0);
		step = _step;
		delay = _delay;
		size = igraph_vcount(&graph);
		isPacemaker.resize(size, false);
		for (int i = 0; i < _np; i++) {
			int tmp;
			file >> tmp;
			isPacemaker[tmp] = true;
		}
		int _nc;
		file >> _nc;
		isContrarian.resize(size, false);
		if (_nc) file >> phaseDif;
		for (int i = 0; i < _nc; i++) {
			int tmp;
			file >> tmp;
			isContrarian[tmp] = true;
		}
		while (file >> _theta >> _omega && theta.size() < size) {
			theta.push_back(vector<double> ());
			theta[theta.size()-1].push_back(_theta);
			freq.push_back(vector<double> ());
			freq[freq.size()-1].push_back(_omega);
			omega.push_back(_omega);
		}
		file.close();
		if (theta.size() != size)
			throw;
		dist2.resize(theta.size());
		for (int i = 0; i < theta.size(); i++) {
			dist2[i].resize(theta.size());
			igraph_vector_t nid;
			igraph_vector_init (&nid, 0); 
			igraph_neighbors(&graph, &nid, i, IGRAPH_IN);				
			int adj_sz = igraph_vector_size(&nid);
			for (int j = 0; j < theta.size(); j++) {
				dist2[i][j] = 1/_sigma;
			}
			igraph_vector_destroy(&nid);
		}
		nivel.resize(size, -1);
		pnivel.resize(size);
	}

	Kakaroto (string fn, string dfn, double _step, int _delay = 0) {
		string gr = fn, conf, line;
		double _theta, _omega, _t0;
		int size, _np, plo;
		ifstream file;
		//leitura dos arquivos
		gr.append(".gr");
		file.open(gr.c_str());
		//leitura do grafo
		getline(file, line);
		Util::readGraph (&graph, line.c_str());
		getline(file, conf);
		cerr << conf << endl;
		file.close();
		//leitura da configuracao
		file.open(conf.c_str());
		file >> _t0 >>_np;
		t.push_back(_t0);
		step = _step;
		delay = _delay;
		size = igraph_vcount(&graph);
		isPacemaker.resize(size, false);
		for (int i = 0; i < _np; i++) {
			int tmp;
			file >> tmp;
			isPacemaker[tmp] = true;
		}
		while (file >> _theta >> _omega && theta.size() < size) {
			theta.push_back(vector<double> ());
			theta[theta.size()-1].push_back(_theta);
			omega.push_back(_omega);
		}
		file.close();
		if (theta.size() != size)
			throw;
		file.open(dfn.c_str());
		dist2.resize(theta.size());
		for (int i = 0; i < theta.size(); i++) {
			dist2[i].resize(theta.size());
			for (int j = 0; j < theta.size(); j++) {
				file >> dist2[i][j];
			}
		}
		file.close();
		nivel.resize(size, -1);
		pnivel.resize(size);
	}
	void connectPacemakersAll ()
	{
		int size = igraph_vcount(&graph);
		for (int i = 0; i < size; ++i)
		{
			if (isPacemaker[i])
			{
				igraph_vector_t nid;
				igraph_vector_init (&nid, 0); 
				igraph_neighbors(&graph, &nid, i, IGRAPH_IN);				
				int adj_sz = igraph_vector_size(&nid);
				std::set<int> N;
				for (int j = 0; j < adj_sz; ++j)
				{
					int next = (int)VECTOR(nid)[j];
					N.insert(next);
				}	
				for (int j = 0; j < size; ++j)
				{
					if (N.count (j) == 0) 
						igraph_add_edge (&graph, i, j);
				}
				igraph_vector_destroy(&nid);
			}
		}
	}

	void faseMediaPorNivel ()
	{
		printf("Niveis\n");
		int size = igraph_vcount(&graph);
		for (int n = 0; n < size; ++n)
		{
			if (pnivel[n].size() == 0) continue;
			double phase = 0.0;
			double sigma = 0.0;
			for (int x = 0; x < pnivel[n].size(); ++x)
			{
				int i = pnivel[n][x];
				phase += theta[i][theta[i].size()-1]/pnivel[n].size();
			}
			for (int x = 0; x < pnivel[n].size(); ++x)
			{
				int i = pnivel[n][x];
				sigma += pow(theta[i][theta[i].size()-1] - phase, 2);
			}
			sigma  = sqrt(sigma/pnivel[n].size());
			printf("%d %.5lf %.5lf\n", n+1, phase, sigma);
		}
	}

	void dumpAll ()
	{
		for (int i = 0; i < theta.size(); ++i)
		{
			stringstream name;
			name << "./dump/" << i << ".theta";
			string st = name.str();
			ofstream dump(st.c_str(), std::ofstream::out);
			for (int j = 0; j < theta[i].size(); ++j)
				dump << 0.0+step*j << " " << theta[i][j] << " " << endl;
			dump.close();			
		}
		for (int i = 0; i < freq.size(); ++i)
		{
			stringstream name;
			name << "./dump/" << i << ".freq";
			string st = name.str();
			ofstream dump(st.c_str(), std::ofstream::out);
			for (int j = 0; j < freq[i].size(); ++j)
				dump << 0.0+step*j << " " << freq[i][j] << " " << endl;
			dump.close();			
		}
	}	

	void niveis ()
	{
		queue<int> q;
		int size = igraph_vcount(&graph);
		for (int i = 0; i < size; ++i) {
			if (isPacemaker[i] || isContrarian[i]) {
				q.push(i);
				nivel[i] = 0;
				pnivel[0].push_back(i);				
			}
		}
		while (!q.empty())
		{
			int curr = q.front(); q.pop();
			igraph_vector_t nid;
			igraph_vector_init (&nid, 0); 
			igraph_neighbors(&graph, &nid, curr, IGRAPH_OUT);				
			int adj_sz = igraph_vector_size(&nid);
			for (int j = 0; j < adj_sz; ++j)
			{
				int next = (int)VECTOR(nid)[j];
				if (nivel[next] != -1) continue;
				nivel[next] = nivel[curr] + 1;
				pnivel[nivel[next]].push_back(next);
				q.push(next);
			}
			igraph_vector_destroy(&nid);
		}
	}

	void draw_niveis ()
	{
		sf::RenderWindow window(sf::VideoMode(800, 600), "My window");
		sf::View view(
			sf::Vector2f(0.0, 0.0), 
			sf::Vector2f(400,300) );
		window.setView(view);

		for (int i = 0; i < theta[0].size() && window.isOpen(); ++i) {
			sf::Event event;
			
			while (window.pollEvent(event))
			{
				// "close requested" event: we close the window
				if (event.type == sf::Event::Closed)
					window.close();
			}
			window.clear(sf::Color::Black);

			int size = igraph_vcount(&graph);	
			vector<double> x(size), y(size);

			double rho = 20.0;
			if (pnivel[0].size() == 1) rho = 0;
			double rinc = 20.0;
			for (int j = 0; j < size; ++j) 
			{
				if (pnivel[j].size() == 0) continue;
				double angle = 0.0;
				double step = 2*M_PI/(double)(pnivel[j].size());
				
				for (int k = 0; k < pnivel[j].size(); ++k) {
					int p = pnivel[j][k];
					double tt = theta[p][i];
					x[p] = rho * cos(angle);
					y[p] = rho * sin(angle);
					angle += step;
				}

				for (int k = 0; k < pnivel[j].size(); k++) 
				{
					int p = pnivel[j][k];
					sf::CircleShape sp(2.0);
					double tt = theta[p][i];
					// sp.setOrigin(2, 2);
					sp.setPosition(x[p]-2, y[p]-2);
					if ( fabs(tt) < 1e-2 
						|| fabs (tt - 2*M_PI) < 1e-2
						|| fabs (tt - M_PI) < 1e-2) sp.setFillColor(sf::Color::Yellow);
					else sp.setFillColor (sf::Color (40, 40, 40) );
					// sp.setFillColor(RGB_from_freq(tt));
					window.draw(sp);
				}
				rho += rinc;
			}

			sf::Vertex A = sf::Vertex( sf::Vector2f(rho+6, rho) );
			sf::Vertex B = sf::Vertex( sf::Vector2f(rho+6, -rho));
			A.color = sf::Color::Blue;
			B.color = sf::Color::Red;
			sf::Vertex line [] = { A, B };
			window.draw(line, 2, sf::Lines);
			
			sf::CircleShape sp(1.0);
			sp.setPosition(rho+6-1.0, 2*rho*(1-R[i])-rho-1.0);
			sp.setFillColor(sf::Color::White);
			window.draw(sp);
			window.display();
			sf::sleep(sf::seconds(0.01));
		}
	}

	double f (int curr, vector<double> k, double coef) {
		double sum = 0;
		int i, size, adj_size, next;
		igraph_vector_t nid;

		size = igraph_vcount(&graph);
		igraph_vector_init (&nid, 0);
		igraph_neighbors(&graph, &nid, curr, IGRAPH_IN);

		adj_size = igraph_vector_size(&nid);

		if (isPacemaker[curr])	return omega[curr];

		for (i = 0; i < adj_size; i++) {
			next = (int)VECTOR(nid)[i];
			if (!isContrarian[curr])
			{
				if (theta[next].size()-1-delay >= 0) {
					sum += 1/dist2[curr][next] * sin ((theta[next][max(0, (int)(theta[next].size()-1-delay))]
						+k[next]*coef)-(theta[curr][theta[curr].size()-1]+k[curr]*coef));
				}	
			}
			else 
			{
				if (theta[next].size()-1-delay >= 0 && !isContrarian[next]) {
					sum += 1/dist2[curr][next] * sin ((theta[next][max(0, (int)(theta[next].size()-1-delay))]
						+k[next]*coef)-(theta[curr][theta[curr].size()-1]+k[curr]*coef) - phaseDif);
				}		
			}
		}
		igraph_vector_destroy(&nid);
		if (stype == 0) return omega[curr]+sum;
		else return adj_size+sum;
	}
	
	double f2 (int curr, vector<double> k, double coef) {
		double sum = 0;
		int i, size, next;
		return theta[curr][theta[curr].size()-1]+k[curr]*coef;	
	}
	double  circle (double ans)
	{
		double aa;
		aa = floor( fabs(ans)/(2*M_PI) );
		if (ans < 0)	ans += aa+2*M_PI;
		else	ans -= aa*2*M_PI;
		return ans;
	}

	void calc (int iter) {
		int it, i, size;
		vector<double> k0, k1, k2, k3, k4;
		double  ans, aa;
		size = igraph_vcount (&graph);
		k0 = vector<double> (size, 0);
		for (it = 1; it <= iter; it++) {
			k1.clear();
			k2.clear();
			k3.clear();
			k4.clear();
			for (i = 0; i < size; i++)
				k1.push_back(f(i, k0, 0));
			for (i = 0; i < size; i++)
				k2.push_back(f(i, k1, step/2));
			for (i = 0; i < size; i++)
				k3.push_back(f(i, k2, step/2));
			for (i = 0; i < size; i++)
				k4.push_back(f(i, k3, step));
			for (i = 0; i < size; i++) {
				ans = theta[i][it-1] + (step/6.0)*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
				while (ans < -M_PI)	ans += 2*M_PI;
				while (ans > M_PI)	ans -= 2*M_PI;
				freq[i].push_back((1.0/6.0)*(k1[i]+2*k2[i]+2*k3[i]+k4[i]));
				theta[i].push_back (ans);
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
			ang.push_back(atan2(r2, r1));
			R1.push_back(r1);
			R2.push_back(r2);
		}
	}

	void writeR (string fn) {
		ofstream file;
		file.open(fn.c_str());
		Util::printRvector(file, R, "R");
		file.close();
	}

	void calcVar (vector<double> & ans) {
		double mean, var, curr;
		int n, t, i;
		ans.resize(theta[0].size());
		for (t = 0; t < theta[0].size(); t++) {
			mean = 0;
			n = 0;
			for (i = 0; i < theta.size(); i++) {
				if (isPacemaker[i]) {
					mean += theta[i][t];
					n++;
				}
			}
			mean /= n;
			var = 0;
			for (i = 0; i < theta.size(); i++) {
				if (!isPacemaker[i]) {
					curr = min (fabs((theta[i][t]-mean)), fabs((mean-theta[i][t]+2*M_PI)));
					var += curr*curr;
				}
			}
			ans[t] = var/(theta.size()-n);
		}
	}

	void calcVarFreq (vector<double> &ans) {

		double mean, var, curr;
		int n, t, i;
		ans.resize(theta[0].size());
		for (t = 0; t < freq[0].size(); t++) {
			mean = 0;
			n = 0;
			for (i = 0; i < freq.size(); i++) {
				mean += freq[i][t];
				n++;
			}
			mean /= n;
			var = 0;
			for (i = 0; i < freq.size(); i++) {
				// curr = dif(theta[i][t+1], theta[i][t]);
				curr = freq[i][t];
				// cout << curr*(1/step) << endl;
				var += (curr-mean)*(curr-mean);
			}
			ans[t] = var/n;
		}
	}

	void draw_graph (void) {
		sf::RenderWindow window(sf::VideoMode(800, 600), "My window");
		sf::View view(
			sf::Vector2f(0.0, 0.0), 
			sf::Vector2f(400,300) );
		window.setView(view);

		for (int i = 0; i < theta[0].size() && window.isOpen(); ++i) {
			sf::Event event;
			
			while (window.pollEvent(event))
			{
				// "close requested" event: we close the window
				if (event.type == sf::Event::Closed)
					window.close();
			}
			window.clear(sf::Color::Black);


			int size = igraph_vcount(&graph);	
			double rho = 120;
			double angle = 0.0;
			double step = 2*M_PI/(double)size;
			
			vector<double> x(size), y(size);
			for (int j = 0; j < size; ++j) {
				double tt = theta[j][i];
				x[j] = rho * cos(tt);
				y[j] = rho * sin(tt);
				// x[j] = rho * cos(angle);
				// y[j] = rho * sin(angle);
				angle += step;
			}

			for (int k = 0; k < size; k++) {
				igraph_vector_t nid;
				igraph_vector_init (&nid, 0);
				igraph_neighbors(&graph, &nid, k, IGRAPH_IN);				
				int adj_sz = igraph_vector_size(&nid);
				for (int j = 0; j < adj_sz; ++j) {
					int next = (int)VECTOR(nid)[j];
					// cout << x[i] << " - " << x[next] <<  endl;
					sf::Vertex A = sf::Vertex( sf::Vector2f(x[k], y[k]) );
					sf::Vertex B = sf::Vertex( sf::Vector2f(x[next], y[next]));
					double ta = theta[k][i]; 
					double tb = theta[k][i]; 
					A.color =  RGB_from_freq(ta); //sf::Color(ta, ta, ta);
					B.color =  RGB_from_freq(tb); //sf::Color(tb, tb, tb);
					sf::Vertex line [] = { A, B };
					window.draw(line, 2, sf::Lines);
				}
				igraph_vector_destroy(&nid);
			} 
			for (int j = 0; j < size; ++j) {
				if (isPacemaker[j] || isContrarian[j]) continue;
				sf::CircleShape sp(2.0);
				double tt = theta[j][i];
				// sp.setOrigin(2, 2);
				sp.setPosition(x[j]-2, y[j]-2);
				// sp.setPosition(0, 0);
				// sp.setFillColor( sf::Color(tt/2*M_PI * 255, tt/2*M_PI * 255, tt/2*M_PI * 255) );
				sp.setFillColor(RGB_from_freq(tt));
				window.draw(sp);
			}

			for (int j = 0; j < size; ++j) {
				if (!isPacemaker[j] && !isContrarian[j]) continue;
				sf::CircleShape sp(3.0, 3);
				double tt = theta[j][i];
				// sp.setOrigin(2, 2);
				sp.setPosition(x[j]-3.0, y[j]-3.0);
				// sp.setPosition(0, 0);
				// sp.setFillColor( sf::Color(tt/2*M_PI * 255, tt/2*M_PI * 255, tt/2*M_PI * 255) );
				sp.setFillColor(RGB_from_freq(tt));
				window.draw(sp);
			}
			
			sf::Vertex A = sf::Vertex( sf::Vector2f(rho+6, rho) );
			sf::Vertex B = sf::Vertex( sf::Vector2f(rho+6, -rho));
			A.color = sf::Color::Blue;
			B.color = sf::Color::Red;
			sf::Vertex line [] = { A, B };
			window.draw(line, 2, sf::Lines);
			

			sf::CircleShape sp(1.0);
			sp.setPosition(rho+6-1.0, 2*rho*(1-R[i])-rho-1.0);
			sp.setFillColor(sf::Color::White);
			window.draw(sp);

			sf::CircleShape med(3.0);
			med.setPosition(R1[i]-3.0, R2[i]-3.0);
			med.setFillColor(sf::Color::White);
			window.draw(med);

			A = sf::Vertex( sf::Vector2f(0, 0) );
			B = sf::Vertex( sf::Vector2f(R1[i], R2[i]));
			A.color = sf::Color::White;
			B.color = sf::Color::White;
			sf::Vertex lin [] = { A, B };
			window.draw(lin, 2, sf::Lines);

			window.display();
			sf::sleep(sf::seconds(0.025));
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
			for (j = 0; j < size; ++j) {
				if (isPacemaker[j])	continue;
				double tt = theta[j][i];
				sf::CircleShape sp(2);
				sp.setPosition(rho*cos(tt), rho*sin(tt)); 
				sp.setFillColor( sf::Color(255*(j+1)/(double)(size+1),
					255*(j+1)/(double)(size+1), 255*(j+1)/(double)(size+1)) );
				//sp.setFillColor (sf::Color(255, 255, 255));
				window.draw(sp);
			}
			for (j = 0; j < size; ++j) {
				if (!isPacemaker[j])	continue;
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
#endif
