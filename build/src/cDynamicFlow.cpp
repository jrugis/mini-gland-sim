/*
 * cDynamicFlow.cpp
 *
 *  Created on: 02/11/2021
 *      Author: cscott
 */

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <tuple>

#include "h5pp/h5pp.h"

#include "cDynamicFlow.hpp"


cDynamicFlow::cDynamicFlow() : loaded(false), dim(0), tstart(0), tend(0) {
}

cDynamicFlow::cDynamicFlow(const std::string& flow_file) {
  load(flow_file);
}

void cDynamicFlow::load(const std::string& flow_file) {
  // load flow tables from HDF5 file
#ifdef DEBUGDYNAMICFLOW
  std::cout << "<DynamicFlow>: loading flow variables from file: " << flow_file << std::endl;
#endif
  h5pp::File ff(flow_file, h5pp::FilePermission::READONLY);
  t_table = ff.readDataset<Eigen::VectorXd>("/t");
  Cl_table = ff.readDataset<Eigen::VectorXd>("/Cl");
  K_table = ff.readDataset<Eigen::VectorXd>("/K");
  Na_table = ff.readDataset<Eigen::VectorXd>("/Na");
  Q_table = ff.readDataset<Eigen::VectorXd>("/Q");
  dim = t_table.size();
  tstart = t_table(0);
  tend = t_table(Eigen::last);

  // check they are all the same size
  

  loaded = true;
}

std::tuple<int, double> cDynamicFlow::get_interp_vals(const double t) const {
  // returns a tuple containing:
  //   - index of the lower bound index bracketing time t
  //   - proportion t is through the bracketing interval
  // t must be within in the range [t_table(0), t_table(last)]

  // find index of time array that is greater than t
  auto itup = std::upper_bound(t_table.begin(), t_table.end(), t);
  int iup = itup - t_table.begin();
  // and lower index
  int idown = iup - 1;
#ifdef DEBUGDYNAMICFLOW
  if (iup == dim) iup--;
  std::cout << "Found bracket for " << t << ": " << idown << " (" << t_table(idown) << ") -> " << iup << " (" << t_table(iup) << ")" << std::endl;
#endif

  // calculate proportion along [a, b) interval
  double prop = (t - t_table(idown)) / (t_table(iup) - t_table(idown));

  return std::tuple<int, double>(idown, prop);
}

double cDynamicFlow::get_Cl(const std::tuple<int, double>& interp_vals) const {
  // lower index and proportion along bracket
  int i = std::get<0>(interp_vals);
  double t = std::get<1>(interp_vals);

  // compute value
  double val;
  if (i + 1 >= dim) {  // check if we are right at the end of the interval
    // return last value
    val = Cl_table(Eigen::last);
#ifdef DEBUGDYNAMICFLOW
    std::cout << "Calc val = " << val << std::endl;
#endif
  }
  else {
    // interpolate
    double a = Cl_table(i);
    double b = Cl_table(i + 1); // TODO: check range
//    double val = std::lerp(a, b, t);  // requires C++20
    val = a + t * (b - a);
#ifdef DEBUGDYNAMICFLOW
    std::cout << "Calc val = " << val << " ( " << a <<  ", " << b << ", " << t << " )" << std::endl;
#endif
  }

  return val;
}

double cDynamicFlow::get_K(const std::tuple<int, double>& interp_vals) const {
  // lower index and proportion along bracket
  int i = std::get<0>(interp_vals);
  double t = std::get<1>(interp_vals);

  // compute value
  double val;
  if (i + 1 >= dim) {  // check if we are right at the end of the interval
    // return last value
    val = K_table(Eigen::last);
#ifdef DEBUGDYNAMICFLOW
    std::cout << "Calc val = " << val << std::endl;
#endif
  }
  else {
    // interpolate
    double a = K_table(i);
    double b = K_table(i + 1); // TODO: check range
//    double val = std::lerp(a, b, t);  // requires C++20
    val = a + t * (b - a);
#ifdef DEBUGDYNAMICFLOW
    std::cout << "Calc val = " << val << "( " << a <<  ", " << b << ", " << t << " )" << std::endl;
#endif
  }

  return val;
}

double cDynamicFlow::get_Na(const std::tuple<int, double>& interp_vals) const {
  // lower index and proportion along bracket
  int i = std::get<0>(interp_vals);
  double t = std::get<1>(interp_vals);

  // compute value
  double val;
  if (i + 1 >= dim) {  // check if we are right at the end of the interval
    // return last value
    val = Na_table(Eigen::last);
#ifdef DEBUGDYNAMICFLOW
    std::cout << "Calc val = " << val << std::endl;
#endif
  }
  else {
    // interpolate
    double a = Na_table(i);
    double b = Na_table(i + 1); // TODO: check range
//    double val = std::lerp(a, b, t);  // requires C++20
    val = a + t * (b - a);
#ifdef DEBUGDYNAMICFLOW
    std::cout << "Calc val = " << val << "( " << a <<  ", " << b << ", " << t << " )" << std::endl;
#endif
  }

  return val;
}

double cDynamicFlow::get_Q(const std::tuple<int, double>& interp_vals) const {
  // lower index and proportion along bracket
  int i = std::get<0>(interp_vals);
  double t = std::get<1>(interp_vals);

  // compute value
  double val;
  if (i + 1 >= dim) {  // check if we are right at the end of the interval
    // return last value
    val = Q_table(Eigen::last);
#ifdef DEBUGDYNAMICFLOW
    std::cout << "Calc val = " << val << std::endl;
#endif
  }
  else {
    // interpolate
    double a = Q_table(i);
    double b = Q_table(i + 1); // TODO: check range
//    double val = std::lerp(a, b, t);  // requires C++20
    val = a + t * (b - a);
#ifdef DEBUGDYNAMICFLOW
    std::cout << "Calc val = " << val << "( " << a <<  ", " << b << ", " << t << " )" << std::endl;
#endif
  }

  return val;
}
