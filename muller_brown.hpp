#ifndef INCLUDE_GUARD_MULLER_BROWN_POTENTIAL_HPP
#define INCLUDE_GUARD_MULLER_BROWN_POTENTIAL_HPP

# include <bits/stdc++.h>
# include <Eigen/Core>

using namespace std;
using namespace Eigen;

double muller_brown_potential(double *x);
void mb_gradient(double *grad, double *x, double step=0.01);
void mb_optimizer(double *x, double alpha=0.1, double step=0.01);


#endif // INCLUDE_GUARD_MULLER_BROWN_POTENTIAL_HPP