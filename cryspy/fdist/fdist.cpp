#include <cmath>
using namespace std;
#include "fdist.hpp"

double dist2(double x1, double y1, double z1, double x2, double y2, double z2) {
        double out_dist = pow(x1-x2,2) + pow(y1-y2,2) + pow(z1-z2,2);
        return out_dist;

}
double dist(double x1, double y1, double z1, double x2, double y2, double z2) {
        double tmp_d2 = dist2(x1, y1, z1, x2, y2, z2);
	double di = sqrt(tmp_d2);
        return di;

}
