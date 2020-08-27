#ifndef INCLUDE_GUARD_MULLER_BROWN_POTENTIAL_HPP
#define INCLUDE_GUARD_MULLER_BROWN_POTENTIAL_HPP

# include <bits/stdc++.h>
# include <Eigen/Core>

using namespace std;
using namespace Eigen;

double muller_brown_potential(double *x);
void mb_gradient(double *grad, double *x, double step=0.01);
void mb_optimizer(double *x, double alpha=0.1, double step=0.01);

double muller_brown_potential(double *x) {
  double energy=0;
  double A[] = {-200,-100,-170,15};
  double a[] = {-1,-1,-6.5,0.7};
  double b[] = {0,0,11,0.6};
  double c[] = {-10,-10,-6.5,0.7};

  double x0[][2] = { {1,0},{0,0.5},{-0.5,1.5},{-1,1}};

  for(int i=0;i<4;i++) {
    double val = 0;
    val += a[i]*(x[0]-x0[i][0])*(x[0]-x0[i][0]); 
    val += b[i]*(x[0]-x0[i][0])*(x[1]-x0[i][1]); 
    val += c[i]*(x[1]-x0[i][1])*(x[1]-x0[i][1]); 
    energy += A[i]*exp(val);
  }

  return energy;
} // dim=2


void mb_gradient(double *grad, double *x, double step) {
  int dim=2;
  double *fwd,*bck;
  fwd = (double*)malloc(sizeof(double)*dim);
  bck = (double*)malloc(sizeof(double)*dim);

  for(int i=0;i<dim;i++) {
    for(int j=0;j<dim;j++) {
      fwd[j] = x[j];
      bck[j] = x[j]; 
    }

    fwd[i] += step;
    bck[i] -= step;
    grad[i] = (muller_brown_potential(fwd)-muller_brown_potential(bck))/(2*step);
  }

  free(fwd);
  free(bck);
}


void mb_optimizer(double *x, double alpha, double step) {
  int dim=2, converge;
  double *grad;
  double threshold = 0.001;
  grad = (double*)malloc(sizeof(double)*dim);

  converge = 0;
  while(converge == 0) {
    converge = 1;
    mb_gradient(grad,x);
    for(int i=0;i<dim;i++) {
      x[i] -= alpha*grad[i];
      if(abs(grad[i]) > threshold) { converge *= 0; }
    }

    // printf("%lf %lf %lf %lf %lf\n", muller_brown_potential(x), x[0], x[1], grad[0], grad[1]);
  }

  free(grad);
}

#endif // INCLUDE_GUARD_MULLER_BROWN_POTENTIAL_HPP