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
	double sigma, step;
	int type;
	string fn = "../networks/";
	if (argc < 3) {
		fprintf (stderr, "Usage: %s <1.type(0-sigma, 1-points)> <2.parameter> <3.step> <4.(optional)draw type(1-per level, 2-moving)> <5.(optional)file name>\n", argv[0]);
		return 1;
	}
	if (argc > 5)	fn.append(string(argv[5]));
	else	fn.append("plo");

	int which = 0;
	if (argc > 4) sscanf (argv[4], "%d", &which);

	sscanf (argv[1], "%d", &type);
	sscanf (argv[3], "%lf", &step);
	if (type)	goku = Kakaroto(fn, argv[2], step);
	else {
		sscanf (argv[2], "%lf", &sigma);
		goku = Kakaroto(fn, sigma, step);
	}
	//goku.draw(fn);
	//goku.connectPacemakersAll();
	goku.calc(10000);

	goku.calcR();
	double R = 0.0;
	for (int i = 200; i < goku.R.size(); ++i)
		R += goku.R[i];
	cerr << "R medio: " << R/(goku.R.size()-200) << endl;
	cerr << "Angulo final: " << goku.ang.back() << endl;
	cout << R/(goku.R.size()-200) << endl; 
	if (which == 1) goku.draw_niveis();
	else if (which == 2) goku.draw_graph();
	goku.writeR("../Resultado/waw.r");
	goku.calcVarFreq (var);
	ofstream plo;
	plo.open ("../Resultado/waw.r", std::ofstream::out | std::ofstream::app);
	Util::printRvector(plo, var, "var");
	plo.close();
	return 0;
}
