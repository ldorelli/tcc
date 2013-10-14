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
	string fn = "";
	if (argc < 3) {
		fprintf (stderr, "Usage: %s <1.sigma> <2.step> <3.(optional)draw type:1-per level, 2-moving> <4.(optional)file name>\n", argv[0]);
		return 1;
	}
	if (argc > 4)	fn.append(string(argv[4]));
	else	fn.append("../networks/plo");

	int which = 0;
	if (argc > 3) sscanf (argv[3], "%d", &which);

	sscanf (argv[1], "%lf", &sigma);
	sscanf (argv[2], "%lf", &step);
	goku = Kakaroto(fn, sigma, step);
	//goku.draw(fn);
	//goku.connectPacemakersAll();
	goku.calc(1000);

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
