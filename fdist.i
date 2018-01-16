/* For swig */
%module fdist

%{

#define SWIG_FILE_WITH_INIT
#include "fidst.hpp"
%}

double dist2(double x1, double y1, double z1, double x2, double y2, double z2);
