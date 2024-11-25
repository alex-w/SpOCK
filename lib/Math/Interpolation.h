//==========================================================================
/*
*    Copyright 2020 Sergio De Florio
*    All rigths reserved
*
*    This file is part of SpOCK
* 
*    SpOCK is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation version 3
* 
*    SpOCK is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
*    GNU General Public License for more details.
* 
*    You should have received a copy of the GNU General Public License
*    along with SpOCK. If not, see <https://www.gnu.org/licenses/>.
*
*/
//==========================================================================

#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

#include <string>
#include <VarTypes.h>
      
#include <Eigen/Core>

using namespace std;
using namespace math;
using namespace Eigen;
 
//------------------------------------------------------------------------------
// Group of functions for INTERPOLATION
//------------------------------------------------------------------------------

// Select an interval of interpolation data given the interpolation vector, the interpolation point and the size of the interval
VectorNi<2> IntplInterval(const vector<double> x_vec, double x_intpl, int intpl_size, bool& valid);
VectorNi<2> IntplInterval(const Ref<const VectorXd>& x_vec, double x_intpl, int intpl_size, bool& valid);
// Perform linear interpolation
double LERP(const Ref<const VectorXd>& x_vec, const Ref<const VectorXd>& y_vec, double x_intpl, bool& valid);
// Perform Lagrange interpolation
double LG_Interpolation(const vector<double> x_vec, const vector<double> y_vec, double x_intpl, bool& valid);
double LG_Interpolation(const Ref<const VectorXd>& x_vec, const Ref<const VectorXd>& y_vec, double x_intpl, bool& valid);
// Perform first derivative of Lagrange interpolation
double LG_Interpolation_FD(const Ref<const VectorXd>& x_vec, const Ref<const VectorXd>& y_vec, double x_intpl, bool& valid);
// Perform cubic Hermite interpolation of state vector (x, y, z, dx/dt, dy/dt, dz/dt)
void CH_Interpolation(const int InterPoints, const int InterpStep, const Ref<const MatrixXd>& timeposvel, Ref<MatrixXd> orbstate_interpolated);
// Perform cubic Hermite interpolation of state vector (x, y, z, dx/dt, dy/dt, dz/dt) at a specific time
void CH_Interpolation(double t_intpl, const Ref<const MatrixXd>& timeposvel, VectorNd<6>& orbstate_interpolated);
// Perform quintic Hermite interpolation of state vector (x, y, z, dx/dt, dy/dt, dz/dt) over a time interval and with a specific interpolation step
void QH_Interpolation(const int InterpStep, const Ref<const MatrixXd>& timeposvelacc, Ref<MatrixXd> orbstate_interpolated);
// Perform quintic Hermite interpolation of state vector and acceleration (x, y, z, dx/dt, dy/dt, dz/dt, dx2/dt2, dy2/dt2, dz2/dt2) at a specific time
void QH_Interpolation(double t_intpl, const Ref<const MatrixXd>& timeposvelacc, VectorNd<9>& orbstate_interpolated);

#endif // INTERPOLATION_H_