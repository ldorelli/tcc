#include "igraph.h"
#include <cmath>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <queue>
#include <set>
#include <SFML/Graphics.hpp>
#include "util.hpp"
#include "kakaroto.hpp"

using namespace std;

int main (int argc, char* argv[]) {
	Kakaroto goku;
	vector<double> theta, omega, var;
	int it = 5000;
	double sigma, step;
	string fn = "";

	if (argc < 3) {
		fprintf (stderr, "Usage: %s <1.sigma> <2.step> <3.(optional)draw type:1-per level, 2-moving> <4. delay> <5. type:0-sigma, 1-degree> <6.(optional)file name>\n"
			 , argv[0]);
		return 1;
	}
	if (argc > 6)	fn.append(string(argv[6]));
	else	fn.append("../networks/plo");

	int which = 0;
	int delay = 0;
	int type = 0;
	int contrarians = 0;
	if (argc > 3) sscanf (argv[3], "%d", &which);
	if (argc > 4) sscanf (argv[4], "%d", &delay);
	if (argc > 5) sscanf (argv[5], "%d", &type);
	sscanf (argv[1], "%lf", &sigma);
	sscanf (argv[2], "%lf", &step);

	goku = Kakaroto(fn, sigma, step, delay, type, contrarians);
	//goku.draw(fn);
	//goku.connectPacemakersAll();
	goku.calc(1000);
	goku.niveis();
	goku.calcR();
	double R = 0.0;
	for (int i = 200; i < goku.R.size(); ++i)
		R += goku.R[i];
	cerr << "R medio: " << R/(goku.R.size()-200) << endl;
	cerr << "Angulo final: " << goku.ang.back() << endl;
	goku.calcVarFreq (var);
	cout << R/(goku.R.size()-200) << " " << var.back() << endl; 
	if (which == 1) goku.draw_niveis();
	else if (which == 2) goku.draw_graph();
	goku.writeR("../Resultado/waw.r");
	// ofstream plo;
	// plo.open ("../Resultado/waw.r", std::ofstream::out | std::ofstream::app);
	// Util::printRvector(plo, var, "var");
	// plo.close();

	// vector<double> freqs;
	// for (int  i = 0; i < goku.freq.size(); i++)	freqs.push_back(goku.freq[i][it-1]);
	// plo.open("../Resultado/freqs.r", std::ofstream::out);
	// Util::printRvector(plo, freqs, "freqs");
	// plo.close();
	// goku.faseMediaPorNivel();
	goku.dumpAll();
	
	return 0;
}
