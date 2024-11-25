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

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include <Constants.h>
#include <VarTypes.h>
#include <Transformations.h>
#include <Interpolation.h>

#include <Eigen/Core>
#include <boost/math/interpolators/cubic_hermite.hpp>
#include <boost/math/interpolators/quintic_hermite.hpp>

using namespace std;
using namespace math;
using namespace Eigen;
using namespace constants;
using namespace mathconst;
using namespace astro;
using namespace boost::math::interpolators;

//------------------------------------------------------------------------------
// VectorNi<2> IntplInterval(const vector<double> x_vec, double x_intpl, int intpl_size, bool valid)
//------------------------------------------------------------------------------
/**
 * Select an interval of interpolation data given the interpolation vector, the interpolation point and the size of the interval
 *
 * @param x_vec         Vector of independent variable values
 * @param x_intpl       Interpolation point
 * @param intpl_size    Number of points to be interpolated ( >= 2 )
 * 
 * @return      Vector containing indices of first and last point of interpolation interval
 * @return      Validity flag
 * 
 */
//------------------------------------------------------------------------------   
VectorNi<2> IntplInterval(const vector<double> x_vec, double x_intpl, int intpl_size, bool& valid)
                        {
                        VectorNi<2> IntplInt = VectorNi<2>::Zero();
                        valid = true;
                        
                        int halfsize = intpl_size/2;
                        
                        int vec_size = x_vec.size();
                          
                        for(int i = 0; i < vec_size - 1; i++)
                          {
                          if( x_intpl == x_vec[i] )
                            {
                            IntplInt(0) = i;
                            IntplInt(1) = i;
                            return(IntplInt);  
                            }
                          if( x_intpl == x_vec[i+1] )
                            {
                            IntplInt(0) = i + 1;
                            IntplInt(1) = i + 1;
                            return(IntplInt);  
                            }
                          if( x_vec[i] < x_intpl && x_intpl < x_vec[i+1] )
                            {
                            // Case in which there are no enough interpolation data before or after the interpolation point
                            if( ( i - (halfsize - 1) ) < 0 || ( (i + 1) + (halfsize - 1) ) > vec_size )
                              {
                              valid = false;
                              return(IntplInt);
                              }
                            
                            // Identify indices of interpolation interval
                            IntplInt(0) = i - (halfsize - 1);  
                            IntplInt(1) = (i + 1) + (halfsize - 1);
                            return(IntplInt);
                            }
                          }
                        
                        // x_intpl is outside the boundaries of the interpolation data
                        valid = false;
                        return(IntplInt);
                        };
VectorNi<2> IntplInterval(const Ref<const VectorXd>& x_vec, double x_intpl, int intpl_size, bool& valid)
                        {
                        VectorNi<2> IntplInt = VectorNi<2>::Zero();
                        valid = true;
                        
                        int halfsize = intpl_size/2;
                        
                        int vec_size = x_vec.size();
                          
                        for(int i = 0; i < vec_size - 1; i++)
                          {
                          if( x_intpl == x_vec(i) )
                            {
                            IntplInt(0) = i;
                            IntplInt(1) = i;
                            return(IntplInt);  
                            }
                          if( x_intpl == x_vec(i+1) )
                            {
                            IntplInt(0) = i + 1;
                            IntplInt(1) = i + 1;
                            return(IntplInt);  
                            }
                          if( x_vec(i) < x_intpl && x_intpl < x_vec(i+1) )
                            {
                            // Case in which there are no enough interpolation data before or after the interpolation point
                            if( ( i - (halfsize - 1) ) < 0 || ( (i + 1) + (halfsize - 1) ) > vec_size )
                              {
                              valid = false;
                              return(IntplInt);
                              }
                            
                            // Identify indices of interpolation interval
                            IntplInt(0) = i - (halfsize - 1);  
                            IntplInt(1) = (i + 1) + (halfsize - 1);
                            return(IntplInt);
                            }
                          }
                        
                        // x_intpl is outside the boundaries of the interpolation data
                        valid = false;
                        return(IntplInt);
                        };
//------------------------------------------------------------------------------
// double LERP(const Ref<const VectorXd>& x_vec, const Ref<const VectorXd>& y_vec, double x_intpl, bool& valid)
//------------------------------------------------------------------------------
/**
 * Perform Lagrange interpolation
 *
 * @param x_vec         Vector of independent variable values
 * @param y_vec         Vector of interpolated function values y = f(x)
 * @param x_intpl       Interpolation point
 * 
 * @return      Interpolated point
 * @return      Validity flag
 * 
 */
//------------------------------------------------------------------------------   
double LERP(const Ref<const VectorXd>& x_vec, const Ref<const VectorXd>& y_vec, double x_intpl, bool& valid)
                        {
                        double y_intpl = 0.0;
                        valid = true;
                        
                        y_intpl = y_vec(0) + ( x_intpl - x_vec(0) )*( y_vec(1) - y_vec(0) )/( x_vec(1) - x_vec(0) );
                        
                        return(y_intpl);  
                        };  
//------------------------------------------------------------------------------
// double LG_Interpolation(const vector<double> x_vec, const vector<double> y_vec, double x_intpl, bool& valid)
//------------------------------------------------------------------------------
/**
 * Perform Lagrange interpolation
 *
 * @param x_vec         Vector of independent variable values
 * @param y_vec         Vector of interpolated function values y = f(x)
 * @param x_intpl       Interpolation point
 * 
 * @return      Interpolated point
 * @return      Validity flag
 * 
 */
//------------------------------------------------------------------------------   
double LG_Interpolation(const vector<double> x_vec, const vector<double> y_vec, double x_intpl, bool& valid)
                        {
                        double y_intpl = 0.0;
                        valid = true;
                        
                        int vec_size = x_vec.size();
                        double L;
                        
                        for(int i = 0; i < vec_size; i++)
                          {
                          L = 1.0;
                          
                          // Build Li(x)
                          for(int j = 0; j < vec_size; j++) if( j != i ) L *= ( x_intpl - x_vec[j] )/( x_vec[i] - x_vec[j] );
                          
                          // Compute interpolated point
                          y_intpl += L*y_vec[i];
                          }
                        
                        return(y_intpl);  
                        }
//------------------------------------------------------------------------------
// double LG_Interpolation(const Ref<const VectorXd>& x_vec, const Ref<const VectorXd>& y_vec, double x_intpl, bool& valid)
//------------------------------------------------------------------------------
/**
 * Perform Lagrange interpolation
 *
 * @param x_vec         Vector of independent variable values
 * @param y_vec         Vector of interpolated function values y = f(x)
 * @param x_intpl       Interpolation point
 * 
 * @return      Interpolated point
 * @return      Validity flag
 * 
 */
//------------------------------------------------------------------------------   
double LG_Interpolation(const Ref<const VectorXd>& x_vec, const Ref<const VectorXd>& y_vec, double x_intpl, bool& valid)
                        {
                        double y_intpl = 0.0;
                        valid = true;
                        
                        int vec_size = x_vec.size();
                        double L;
                        
                        for(int i = 0; i < vec_size; i++)
                          {
                          L = 1.0;
                          
                          // Build Li(x)
                          for(int j = 0; j < vec_size; j++) if(j != i) L *= ( x_intpl - x_vec(j) )/( x_vec(i) - x_vec(j) );
                          
                          // Compute interpolated point
                          y_intpl += L*y_vec(i);
                          }
                        
                        return(y_intpl);  
                        }
//------------------------------------------------------------------------------
// double LG_Interpolation_FD(const Ref<const VectorXd>& x_vec, const Ref<const VectorXd>& y_vec, double x_intpl, bool& valid)
//------------------------------------------------------------------------------
/**
 * Perform Lagrange interpolation
 *
 * @param x_vec         Vector of independent variable values
 * @param y_vec         Vector of interpolated function values y = f(x)
 * @param x_intpl       Interpolation point
 * 
 * @return      Interpolated point
 * @return      Validity flag
 * 
 */
//------------------------------------------------------------------------------   
double LG_Interpolation_FD(const Ref<const VectorXd>& x_vec, const Ref<const VectorXd>& y_vec, double x_intpl, bool& valid)
                        {
                        double y_intpl_FD = 0.0;
                        valid = true;
                        
                        int vec_size = x_vec.size();
                        double L, L_FD;
                        
                        for(int i = 0; i < vec_size; i++)
                          {
                          L_FD = 0.0;
                          
                          // Build L_FDi(x)
                          for(int k = 0; k < vec_size; k++)
                              {
                              L = 1.0;
                              
                              for(int j = 0; j < vec_size; j++) if(j != i) L *= ( x_intpl - x_vec(j) )/( x_vec(i) - x_vec(j) );
                              
                              if(k != i) L_FD += L/(x_intpl - x_vec(k)); // x_intpl == x_vec(k) in case x_intpl coincide with one of the points to be interpolated
                              }
                          
                          // Compute interpolated point
                          y_intpl_FD += L_FD*y_vec(i);
                          }
                        
                        return(y_intpl_FD);  
                        }                        
//------------------------------------------------------------------------------
// void QH_Interpolation(const int InterpStep, const Ref<const MatrixXd>& timeposvelacc, Ref<MatrixXd> orbstate_interpolated)
//------------------------------------------------------------------------------
/**
 * Perform quintic Hermite interpolation of state vector (x, y, z, dx/dt, dy/dt, dz/dt) over a time interval and with a specific interpolation step
 *
 * @param InterpStep            Step of interpolated points                 
 * @param timeposvelacc         Matrix of points used for the computation of the interpolation polynomial
 *                              (each row contains time, x, y, z, vx, vy, vz, ax, ay, az)
 * 
 * @return Matrix of interpolated points (each row contains time, x, y, z, vx, vy, vz)
 */
//------------------------------------------------------------------------------   
void QH_Interpolation(const int InterpStep,
                      const Ref<const MatrixXd>& timeposvelacc,
                      Ref<MatrixXd> orbstate_interpolated)
                      {
                      const int InterPoints = timeposvelacc.rows();
                      double t_intpl;
                      int ind = 1;
                      
                      vector<double> v_null(InterPoints);
                      for(int i = 0; i < InterPoints; i++) v_null[i] = 0.0;
                      
                      vector<double> time_interpol(InterPoints);
                      
                      vector<double> x_interpol(InterPoints);
                      vector<double> dxdt_interpol(InterPoints);
                      vector<double> dx2dt2_interpol(InterPoints);
                      
                      int lastrow = orbstate_interpolated.rows() - 1;
                      
                      orbstate_interpolated.row(0) = timeposvelacc.row(0).segment(0,7);
                      
                      for(int i = 0; i < InterPoints; i++)
                          {
                          time_interpol[i] = timeposvelacc(i,0);
                          x_interpol[i] = timeposvelacc(i,1);
                          dxdt_interpol[i] = timeposvelacc(i,4);
                          dx2dt2_interpol[i] = timeposvelacc(i,7);
                          }
                      
                      auto spline_x = quintic_hermite(std::move(time_interpol), std::move(x_interpol), std::move(dxdt_interpol), std::move(dx2dt2_interpol));
                      
                      time_interpol.clear();   time_interpol = v_null;
                      x_interpol.clear();      x_interpol = v_null;
                      dxdt_interpol.clear();   dxdt_interpol = v_null;
                      dx2dt2_interpol.clear(); dx2dt2_interpol = v_null;
                      
                      for(int i = 0; i < InterPoints; i++)
                          {
                          time_interpol[i] = timeposvelacc(i,0);
                          x_interpol[i] = timeposvelacc(i,2);
                          dxdt_interpol[i] = timeposvelacc(i,5);
                          dx2dt2_interpol[i] = timeposvelacc(i,8);
                          }
                      
                      auto spline_y = quintic_hermite(std::move(time_interpol), std::move(x_interpol), std::move(dxdt_interpol), std::move(dx2dt2_interpol));
                      
                      time_interpol.clear();   time_interpol = v_null;
                      x_interpol.clear();      x_interpol = v_null;
                      dxdt_interpol.clear();   dxdt_interpol = v_null;
                      dx2dt2_interpol.clear(); dx2dt2_interpol = v_null;
                      
                      for(int i = 0; i < InterPoints; i++)
                          {
                          time_interpol[i] = timeposvelacc(i,0);
                          x_interpol[i] = timeposvelacc(i,3);
                          dxdt_interpol[i] = timeposvelacc(i,6);
                          dx2dt2_interpol[i] = timeposvelacc(i,9);
                          }
                      
                      auto spline_z = quintic_hermite(std::move(time_interpol), std::move(x_interpol), std::move(dxdt_interpol), std::move(dx2dt2_interpol));
                      
                      time_interpol.clear();   time_interpol = v_null;
                      x_interpol.clear();      x_interpol = v_null;
                      dxdt_interpol.clear();   dxdt_interpol = v_null;
                      dx2dt2_interpol.clear(); dx2dt2_interpol = v_null;
                      
                      
                      t_intpl = orbstate_interpolated(0,0) + InterpStep;
                      
                      for(int i = 1; i < lastrow; i++)
                          {
                          if( t_intpl == timeposvelacc(ind,0) )
                              {
                               orbstate_interpolated.row(i) = timeposvelacc.row(ind).segment(0,7);
                               
                               t_intpl += InterpStep;
                               ind++;
                               //continue;
                              }
                          else
                              {
                              orbstate_interpolated(i,0) = t_intpl;
                              orbstate_interpolated(i,1) = spline_x(t_intpl);
                              orbstate_interpolated(i,2) = spline_y(t_intpl);
                              orbstate_interpolated(i,3) = spline_z(t_intpl);
                              orbstate_interpolated(i,4) = spline_x.prime(t_intpl);
                              orbstate_interpolated(i,5) = spline_y.prime(t_intpl);
                              orbstate_interpolated(i,6) = spline_z.prime(t_intpl);
                              
                              t_intpl += InterpStep;
                              }
                          }
                          
                      orbstate_interpolated.row(lastrow) = timeposvelacc.row(InterPoints-1).segment(0,7);
                      };
//------------------------------------------------------------------------------
// void QH_Interpolation(double t_intpl, const Ref<const MatrixXd>& timeposvelacc, VectorNd<9>& orbstate_interpolated)
//------------------------------------------------------------------------------
/**
 * Perform quintic Hermite interpolation of state vector and acceleration (x, y, z, dx/dt, dy/dt, dz/dt, dx2/dt2, dy2/dt2, dz2/dt2) at a specific time
 *
 * @param t_intpl               Interpolation time
 * @param timeposvelacc         Matrix of points used for the computation of the interpolation polynomial
 *                              (each row contains time, x, y, z, vx, vy, vz, ax, ay, az)
 * 
 * @return Interpolated state vector and acceleration (x, y, z, dx/dt, dy/dt, dz/dt, dx2/dt2, dy2/dt2, dz2/dt2)
 */
//------------------------------------------------------------------------------   
void QH_Interpolation(double t_intpl,
                      const Ref<const MatrixXd>& timeposvelacc,
                      VectorNd<9>& orbstate_interpolated)
                      {
                      const int InterPoints = timeposvelacc.rows();
                      
                      for(int i = 0; i < InterPoints; i++)
                          {
                          if(t_intpl == timeposvelacc(i,0)) // No need of interpolation
                              {
                              orbstate_interpolated = timeposvelacc.row(i).segment(1,9);
                              return;
                              }
                          }
                      
                      vector<double> v_null(InterPoints);
                      for(int i = 0; i < InterPoints; i++) v_null[i] = 0.0;
                      
                      vector<double> time_interpol(InterPoints);
                      vector<double> x_interpol(InterPoints);
                      vector<double> dxdt_interpol(InterPoints);
                      vector<double> dx2dt2_interpol(InterPoints);
                      
                      for(int i = 0; i < InterPoints; i++)
                          {
                          time_interpol[i] = timeposvelacc(i,0);
                          x_interpol[i] = timeposvelacc(i,1);
                          dxdt_interpol[i] = timeposvelacc(i,4);
                          dx2dt2_interpol[i] = timeposvelacc(i,7);
                          }
                      
                      auto spline_x = quintic_hermite(std::move(time_interpol), std::move(x_interpol), std::move(dxdt_interpol), std::move(dx2dt2_interpol));
                      
                      time_interpol.clear();   time_interpol = v_null;
                      x_interpol.clear();      x_interpol = v_null;
                      dxdt_interpol.clear();   dxdt_interpol = v_null;
                      dx2dt2_interpol.clear(); dx2dt2_interpol = v_null;
                      
                      for(int i = 0; i < InterPoints; i++)
                          {
                          time_interpol[i] = timeposvelacc(i,0);
                          x_interpol[i] = timeposvelacc(i,2);
                          dxdt_interpol[i] = timeposvelacc(i,5);
                          dx2dt2_interpol[i] = timeposvelacc(i,8);
                          }
                      
                      auto spline_y = quintic_hermite(std::move(time_interpol), std::move(x_interpol), std::move(dxdt_interpol), std::move(dx2dt2_interpol));
                      
                      time_interpol.clear();   time_interpol = v_null;
                      x_interpol.clear();      x_interpol = v_null;
                      dxdt_interpol.clear();   dxdt_interpol = v_null;
                      dx2dt2_interpol.clear(); dx2dt2_interpol = v_null;
                      
                      for(int i = 0; i < InterPoints; i++)
                          {
                          time_interpol[i] = timeposvelacc(i,0);
                          x_interpol[i] = timeposvelacc(i,3);
                          dxdt_interpol[i] = timeposvelacc(i,6);
                          dx2dt2_interpol[i] = timeposvelacc(i,9);
                          }
                      
                      auto spline_z = quintic_hermite(std::move(time_interpol), std::move(x_interpol), std::move(dxdt_interpol), std::move(dx2dt2_interpol));
                      
                      time_interpol.clear();   time_interpol = v_null;
                      x_interpol.clear();      x_interpol = v_null;
                      dxdt_interpol.clear();   dxdt_interpol = v_null;
                      dx2dt2_interpol.clear(); dx2dt2_interpol = v_null;
                      
                      // Interpolation                                
                      orbstate_interpolated(0) = spline_x(t_intpl);
                      orbstate_interpolated(1) = spline_y(t_intpl);
                      orbstate_interpolated(2) = spline_z(t_intpl);
                      orbstate_interpolated(3) = spline_x.prime(t_intpl);
                      orbstate_interpolated(4) = spline_y.prime(t_intpl);
                      orbstate_interpolated(5) = spline_z.prime(t_intpl);
                      orbstate_interpolated(6) = spline_x.double_prime(t_intpl);
                      orbstate_interpolated(7) = spline_y.double_prime(t_intpl);
                      orbstate_interpolated(8) = spline_z.double_prime(t_intpl);
                      };
//------------------------------------------------------------------------------
// void CH_Interpolation(const int InterPoints, const int InterpStep, const Ref<const MatrixXd>& timeposvel, Ref<MatrixXd> orbstate_interpolated)
//------------------------------------------------------------------------------
/**
 * Perform cubic Hermite interpolation of state vector (x, y, z, dx/dt, dy/dt, dz/dt)
 *
 * @param InterPoints           Number of points used for the computation of the interpolation polynomial
 * @param InterpStep            Step of interpolated points                 
 * @param timeposvel            Matrix of points used for the computation of the interpolation polynomial
 *                              (each row contains time, x, y, z, vx, vy, vz)
 * 
 * @return Matrix of interpolated points (each row contains time, x, y, z, vx, vy, vz)  
 */
//------------------------------------------------------------------------------   
void CH_Interpolation(const int InterPoints,
                      const int InterpStep,
                      const Ref<const MatrixXd>& timeposvel,
                      Ref<MatrixXd> orbstate_interpolated)
                      {
                      double t_intpl;
                      int ind = 1;
                      
                      vector<double> v_null(InterPoints);
                      for(int i = 0; i < InterPoints; i++) v_null[i] = 0.0;
                      
                      vector<double> time_interpol(InterPoints);
                      
                      vector<double> x_interpol(InterPoints);
                      vector<double> dxdt_interpol(InterPoints);
                      
                      int lastrow = orbstate_interpolated.rows() - 1;
                      
                      orbstate_interpolated.row(0) = timeposvel.row(0).segment(0,7);
                      
                      //cout << orbstate_interpolated(0,0) << "   " << orbstate_interpolated(0,1) << "   " << orbstate_interpolated(0,2) << "   " << orbstate_interpolated(0,3) << "   " << orbstate_interpolated(0,4) << "   " << orbstate_interpolated(0,5) << "   " << orbstate_interpolated(0,6) << endl;
                      
                      for(int i = 0; i < InterPoints; i++)
                          {
                          time_interpol[i] = timeposvel(i,0);
                          x_interpol[i] = timeposvel(i,1);
                          dxdt_interpol[i] = timeposvel(i,4);
                          }
                      
                      auto spline_x = cubic_hermite(std::move(time_interpol), std::move(x_interpol), std::move(dxdt_interpol));
                      
                      time_interpol.clear();   time_interpol = v_null;
                      x_interpol.clear();      x_interpol = v_null;
                      dxdt_interpol.clear();   dxdt_interpol = v_null;
                      
                      for(int i = 0; i < InterPoints; i++)
                          {
                          time_interpol[i] = timeposvel(i,0);
                          x_interpol[i] = timeposvel(i,2);
                          dxdt_interpol[i] = timeposvel(i,5);
                          }
                      
                      auto spline_y = cubic_hermite(std::move(time_interpol), std::move(x_interpol), std::move(dxdt_interpol));
                      
                      time_interpol.clear();   time_interpol = v_null;
                      x_interpol.clear();      x_interpol = v_null;
                      dxdt_interpol.clear();   dxdt_interpol = v_null;
                      
                      for(int i = 0; i < InterPoints; i++)
                          {
                          time_interpol[i] = timeposvel(i,0);
                          x_interpol[i] = timeposvel(i,3);
                          dxdt_interpol[i] = timeposvel(i,6);
                          }
                      
                      auto spline_z = cubic_hermite(std::move(time_interpol), std::move(x_interpol), std::move(dxdt_interpol));
                      
                      time_interpol.clear();   time_interpol = v_null;
                      x_interpol.clear();      x_interpol = v_null;
                      dxdt_interpol.clear();   dxdt_interpol = v_null;
                      
                      
                      t_intpl = orbstate_interpolated(0,0) + InterpStep;
                      
                      for(int i = 1; i < lastrow; i++)
                          {
                          //cout << "t_intpl = " << t_intpl << endl;
                          //cout << "timeposvel(ind,0) = " << timeposvel(ind,0) << endl;
                          if( t_intpl == timeposvel(ind,0) )
                              {
                               orbstate_interpolated.row(i) = timeposvel.row(ind).segment(0,7);
                               
                               t_intpl += InterpStep;
                               ind++;
                               //continue;
                              }
                          else
                              {
                              orbstate_interpolated(i,0) = t_intpl;
                              
                              orbstate_interpolated(i,1) = spline_x(t_intpl);
                              orbstate_interpolated(i,4) = spline_x.prime(t_intpl);
                              
                              orbstate_interpolated(i,2) = spline_y(t_intpl);
                              orbstate_interpolated(i,5) = spline_y.prime(t_intpl);
                              
                              orbstate_interpolated(i,3) = spline_z(t_intpl);
                              orbstate_interpolated(i,6) = spline_z.prime(t_intpl);
                              
                              //t_intpl++;
                              t_intpl += InterpStep;
                              }
                          
                          //cout << orbstate_interpolated(i,0) << "   " << orbstate_interpolated(i,1) << "   " << orbstate_interpolated(i,2) << "   " << orbstate_interpolated(i,3) << "   " << orbstate_interpolated(i,4) << "   " << orbstate_interpolated(i,5) << "   " << orbstate_interpolated(i,6) << endl;
                          }
                          
                      orbstate_interpolated.row(lastrow) = timeposvel.row(InterPoints-1).segment(0,7);
                          
                      //cout << orbstate_interpolated(lastrow,0) << "   " << orbstate_interpolated(lastrow,1) << "   " << orbstate_interpolated(lastrow,2) << "   " << orbstate_interpolated(lastrow,3) << "   " << orbstate_interpolated(lastrow,4) << "   " << orbstate_interpolated(lastrow,5) << "   " << orbstate_interpolated(lastrow,6) << endl;   
                      };
//------------------------------------------------------------------------------
// void CH_Interpolation(double t_intpl, const Ref<const MatrixXd>& timeposvel, VectorNd<7>& orbstate_interpolated)
//------------------------------------------------------------------------------
/**
 * Perform cubic Hermite interpolation of state vector (x, y, z, dx/dt, dy/dt, dz/dt) at a specific time
 *
 * @param t_intpl               Interpolation time
 * @param timeposvel            Matrix of points used for the computation of the interpolation polynomial
 *                              (each row contains time, x, y, z, vx, vy, vz)
 * 
 * @return Interpolated state vector (x, y, z, vx, vy, vz)
 */
//------------------------------------------------------------------------------   
void CH_Interpolation(double t_intpl,
                      const Ref<const MatrixXd>& timeposvel,
                      VectorNd<6>& orbstate_interpolated)
                      {
                      const int InterPoints = timeposvel.rows();
                      
                      for(int i = 0; i < InterPoints; i++)
                          {
                          if(t_intpl == timeposvel(i,0)) // No need of interpolation
                              {
                              orbstate_interpolated = timeposvel.row(i).segment(1,6);
                              return;
                              }
                          }
                      
                      vector<double> v_null(InterPoints);
                      for(int i = 0; i < InterPoints; i++) v_null[i] = 0.0;
                      
                      vector<double> time_interpol(InterPoints);
                      vector<double> x_interpol(InterPoints);
                      vector<double> dxdt_interpol(InterPoints);
                      
                      for(int i = 0; i < InterPoints; i++)
                          {
                          time_interpol[i] = timeposvel(i,0);
                          x_interpol[i] = timeposvel(i,1);
                          dxdt_interpol[i] = timeposvel(i,4);
                          }
                      
                      auto spline_x = cubic_hermite(std::move(time_interpol), std::move(x_interpol), std::move(dxdt_interpol));
                      
                      time_interpol.clear();   time_interpol = v_null;
                      x_interpol.clear();      x_interpol = v_null;
                      dxdt_interpol.clear();   dxdt_interpol = v_null;
                      
                      for(int i = 0; i < InterPoints; i++)
                          {
                          time_interpol[i] = timeposvel(i,0);
                          x_interpol[i] = timeposvel(i,2);
                          dxdt_interpol[i] = timeposvel(i,5);
                          }
                      
                      auto spline_y = cubic_hermite(std::move(time_interpol), std::move(x_interpol), std::move(dxdt_interpol));
                      
                      time_interpol.clear();   time_interpol = v_null;
                      x_interpol.clear();      x_interpol = v_null;
                      dxdt_interpol.clear();   dxdt_interpol = v_null;
                      
                      for(int i = 0; i < InterPoints; i++)
                          {
                          time_interpol[i] = timeposvel(i,0);
                          x_interpol[i] = timeposvel(i,3);
                          dxdt_interpol[i] = timeposvel(i,6);
                          }
                      
                      auto spline_z = cubic_hermite(std::move(time_interpol), std::move(x_interpol), std::move(dxdt_interpol));
                      
                      time_interpol.clear();   time_interpol = v_null;
                      x_interpol.clear();      x_interpol = v_null;
                      dxdt_interpol.clear();   dxdt_interpol = v_null;
                      
                      // Interpolation                                
                      orbstate_interpolated(0) = spline_x(t_intpl);
                      orbstate_interpolated(1) = spline_y(t_intpl);
                      orbstate_interpolated(2) = spline_z(t_intpl);
                      orbstate_interpolated(3) = spline_x.prime(t_intpl);
                      orbstate_interpolated(4) = spline_y.prime(t_intpl);
                      orbstate_interpolated(5) = spline_z.prime(t_intpl);
                      };
