/*
 * cDynamicFlow.hpp
 *
 *  Created on: 02/11/2021
 *      Author: cscott
 */

#ifndef CDYNAMICFLOW_H_
#define CDYNAMICFLOW_H_

#include <string>
#include <tuple>
#include <Eigen/Core>
#include <Eigen/Dense>

class cDynamicFlow {
  public:
    cDynamicFlow();
    cDynamicFlow(const std::string& flow_file);
    void load(const std::string& flow_file);
    std::tuple<int, double> get_interp_vals(const double t) const;
    double get_Cl(const std::tuple<int, double>& interp_vals) const;
    double get_K(const std::tuple<int, double>& interp_vals) const;
    double get_Na(const std::tuple<int, double>& interp_vals) const;
    double get_Q(const std::tuple<int, double>& interp_vals) const;
    double get_tstart() const { return tstart; }
    double get_tend() const { return tend; }

  private:
    bool loaded;
    int dim;  // number of points in the table
    double tstart, tend;  // start and end times of the input
    Eigen::VectorXd t_table;
    Eigen::VectorXd Cl_table;
    Eigen::VectorXd K_table;
    Eigen::VectorXd Na_table;
    Eigen::VectorXd Q_table;
};

#endif /* CDYNAMICFLOW_H_ */
