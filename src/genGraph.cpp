#include "igraph.h"
#include "util.hpp"

int main (int argc, char* argv[]) {
	if (argc < 5) {
		printf ("Usage: %s <1.population> <2.type: 1-BA 2-ER 3-WS 4-NLBA 5-SF2ER> <3.parameters> <4.output file>\n", argv[0]);
		return 1;
	}
	srand(time(NULL));
	Util::genGraph(atoi(argv[1]), atoi(argv[2]), argv[3], argv[4]);
	return 0;
}
