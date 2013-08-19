#include <iostream>
#include <cstdlib>
#include <cmath>
// Standard
double normal(double average, double stdev) {		
	double r = sqrt(-2.0 * log(random()/(RAND_MAX * 1.0)));
	double theta = 2.0 * M_PI * random() / (RAND_MAX + 1.0);
	double val = average + stdev * r * cos(theta);
	return val;
}
