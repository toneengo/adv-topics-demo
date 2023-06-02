#ifndef COMP_CPP
#define COMP_CPP

#include <math.h>

class Comp {
    public:
        static bool dcomp(double a, double b, double tol=1e-6) { 
            return (abs(a - b) < tol ? true : false); 
        }
};

#endif