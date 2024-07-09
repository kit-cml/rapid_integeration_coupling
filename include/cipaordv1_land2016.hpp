#ifndef CIPAORDV1_LAND2016_HPP
#define CIPAORDV1_LAND2016_HPP

#include "Cellmodel.h"
#include "enum_cipaordv1_land2016.hpp"

#define EULER

class cipaordv1_land2016 : public Cellmodel
{
public:
  cipaordv1_land2016();
  ~cipaordv1_land2016();
  void initConsts ();
  void initConsts(double type);
  // void initConsts(double type, double conc, const double *hill, const double *herg );
  // void initConsts(bool is_skinned, bool BETA);
  void initConsts(double type, double conc, const double *hill, const double *herg, bool is_skinned, bool BETA );
  void computeRates( double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC );
  // void computeRates( double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC, double land_trpn );
  void solveEuler(double dt);
  // void solveAnalytical( int forward_euler_only, double dt, double *CONSTANTS, double *RATES, double* STATES, double* ALGEBRAIC );
  double tryDt( double dt,
               double TIME,
               double *CONSTANTS,
               double *RATES,
               double* STATES,
               double* ALGEBRAIC );
  // void gaussElimination(double *A, double *b, double *x, int N);
//  void gaussSeidel(double **a, double *b, double *x, int n, int maxIterations, double tolerance);
  void solveRK4(double TIME,double dt);
  // void mat_vec_multiply(double** a, double* x, double *result, int n);
  // void solve_rk_hERG(double** a, double* y, double dt, int n );
  double set_time_step(double TIME,
                       double time_point,
                       double min_time_step,
                       double max_time_step,
                       double min_dV,
                       double max_dV,
                       double* CONSTANTS,
                       double* RATES,
                       double* STATES,
                       double* ALGEBRAIC);
private:
  void ___applyDrugEffect(double conc, const double *hill);
  void ___applyHERGBinding(double conc, const double *herg);
  void ___initConsts(double type);
  void ___initConsts(bool is_skinned,
                     bool BETA);
};


#endif

