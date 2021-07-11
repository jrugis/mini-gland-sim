/*
 * utils.hpp
 *
 *  Created on: 21/04/2021
 *      Author: jrugis
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "inih/cpp/INIReader.h"
#include "global_defs.hpp"

namespace utils {
  void calc_tri_centers(MatrixN3d& centers, const MatrixN3d& vertices, const MatrixN3i& triangles);
  void calc_tri_areas(MatrixN1d& areas, const MatrixN3d& vertices, const MatrixN3i& triangles);
  void fatal_error(const std::string msg, std::ofstream& out);
  void get_parameters(const std::string file_name, std::unordered_map<std::string, double>& p, std::ofstream& out);
  void assert_parameter_set(const double value, const std::string& error_msg, std::ofstream& out);
  void assert_parameter_set(const int value, const std::string& error_msg, std::ofstream& out);
  void assert_parameter_set(const std::string& value, const std::string& error_msg, std::ofstream& out);
  double get_parameter_real(INIReader* ini, const std::string& section, const std::string& name, std::ofstream& out);
  int get_parameter_integer(INIReader* ini, const std::string& section, const std::string& name, std::ofstream& out);
  std::string get_parameter_string(INIReader* ini, const std::string& section, const std::string& name, std::ofstream& out);
  double get_distance(const Vector3d& p, const Vector3d& v, const Vector3d& w);
  void mesh_open(std::ifstream& mesh_file, std::ofstream& out);
  int mesh_get_count(std::ifstream& file, std::string tag);
  void mesh_end_header(std::ifstream& file);
  void mesh_skip_lines(std::ifstream& file, int count);
  std::vector<std::string> mesh_get_tokens(std::ifstream& file);
  void save_matrix(const std::string file_name, int bytes, char* data);
} // namespace utils

#endif /* UTILS_H_ */
