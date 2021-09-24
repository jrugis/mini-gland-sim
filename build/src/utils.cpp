/*
 * utils.cpp
 *
 *	Created on: 21/04/2021
 *	Author: jrugis
 */

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/tokenizer.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <cmath>
#include <fstream>
#include <iostream>

#include "global_defs.hpp"
#include "utils.hpp"

// get the center points of triangles
void utils::calc_tri_centers(MatrixN3d& centers, const MatrixN3d& vertices, const MatrixN3i& triangles)
{
  centers.resize(triangles.rows(), Eigen::NoChange);
  for (int n = 0; n < triangles.rows(); n++) {
    Vector3d v1 = Vector3d(vertices.row((triangles)(n, 0)));
    Vector3d v2 = Vector3d(vertices.row((triangles)(n, 1)));
    Vector3d v3 = Vector3d(vertices.row((triangles)(n, 2)));
    centers.row(n) = (v1 + v2 + v3) / 3.0;
  }
}

// calculate the area of each triangle
void utils::calc_tri_areas(MatrixN1d& areas, const MatrixN3d& vertices, const MatrixN3i& triangles)
{
  areas.resize(triangles.rows(), Eigen::NoChange);
  for (int n = 0; n < triangles.rows(); n++) {
    Vector3d v1 = Vector3d(vertices.row((triangles)(n, 0)));
    Vector3d v2 = Vector3d(vertices.row((triangles)(n, 1)));
    Vector3d v3 = Vector3d(vertices.row((triangles)(n, 2)));
    // Eigen::Matrix<double, 1, 3> side1 = vert.block<1, 3>(0, 0) - vert.block<1, 3>(1, 0);
    // Eigen::Matrix<double, 1, 3> side2 = vert.block<1, 3>(0, 0) - vert.block<1, 3>(2, 0);
    Vector3d side1 = v1 - v2;
    Vector3d side2 = v1 - v3;
    areas(n) = 0.5 * (side1.cross(side2)).norm();
  }
}

// NOTE: outputs error message to stderr and the process "out" file
void utils::fatal_error(const std::string msg, std::ofstream& out)
{
  std::string m = "ERROR: " + msg;
  out << m << std::endl;
  out.close();
  std::cerr << m << std::endl;
  exit(1);
}

// utility functions for getting parameters and at the same time checking they were actually set in the ini file
double utils::get_parameter_real(INIReader* ini, const std::string& section, const std::string& name, std::ofstream& out) {
  double value = ini->GetReal(section, name, PARAMETER_DEFAULT_REAL);
  assert_parameter_set(value, std::string("Parameter not set in input file: ") + section + std::string("::") + name, out);
  return value;
}

int utils::get_parameter_integer(INIReader* ini, const std::string& section, const std::string& name, std::ofstream& out) {
  int value = ini->GetInteger(section, name, PARAMETER_DEFAULT_INTEGER);
  assert_parameter_set(value, std::string("Parameter not set in input file: ") + section + std::string("::") + name, out);
  return value;
}

std::string utils::get_parameter_string(INIReader* ini, const std::string& section, const std::string& name, std::ofstream& out) {
  std::string value = ini->Get(section, name, PARAMETER_DEFAULT_STRING);
  assert_parameter_set(value, std::string("Parameter not set in input file: ") + section + std::string("::") + name, out);
  return value;
}

void utils::assert_parameter_set(const double value, const std::string& error_msg, std::ofstream& out) {
  if (value == PARAMETER_DEFAULT_REAL) {
    fatal_error(error_msg, out);
  }
}

void utils::assert_parameter_set(const int value, const std::string& error_msg, std::ofstream& out) {
  if (value == PARAMETER_DEFAULT_INTEGER) {
    fatal_error(error_msg, out);
  }
}

void utils::assert_parameter_set(const std::string& value, const std::string& error_msg, std::ofstream& out) {
  if (value == PARAMETER_DEFAULT_STRING) {
    fatal_error(error_msg, out);
  }
}

double utils::get_distance(const Vector3d& p, const Vector3d& v, const Vector3d& w)
{
  // Return minimum distance between line segment vw and point p
  double l2 = (w - v).squaredNorm();      // |w-v|^2		avoid a sqrt
  if (l2 == 0.0) return ((v - p).norm()); // v == w case, return distance(p, v)
  // Consider the line extending the segment, parameterized as v + t (w - v).
  // Find projection of point p onto the line. It falls where t = [(p-v) . (w-v)] / |w-v|^2
  // Clamp t from [0,1] to handle points outside the segment vw.
  double t = std::max(0.0, std::min(1.0, (p - v).dot(w - v) / l2)); // max(0, min(1, dot(p - v, w - v) / l2));
  const Vector3d projection = v + (t * (w - v));                    // Projection falls on the segment
  return ((projection - p).norm());                                 // return distance(p, projection)
}

void utils::mesh_open(std::ifstream& mesh_file, std::string mesh_file_name, std::ofstream& out)
{
  mesh_file.open(mesh_file_name); // open the mesh file
  if (not mesh_file.is_open()) fatal_error("mesh file " + std::string(MESH_FILE_NAME) + " could not be opened", out);
  return;
}

int utils::mesh_get_count(std::ifstream& mesh_file, std::string tag)
{
  std::string line;                // file line buffer
  std::vector<std::string> tokens; // tokenized line
  while (getline(mesh_file, line)) {
    boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
    if (tokens[1] == tag) break;
  }
  return (std::stoi(tokens[2]));
}

void utils::mesh_end_header(std::ifstream& mesh_file)
{
  std::string line;                // file line buffer
  std::vector<std::string> tokens; // tokenized line
  while (getline(mesh_file, line)) {
    boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
    if (tokens[0] == std::string("end_header")) break;
  }
  return;
}

void utils::mesh_skip_lines(std::ifstream& mesh_file, int count)
{
  std::string line; // file line buffer
  for (int i = 0; i < count; i++) getline(mesh_file, line);
  return;
}

std::vector<std::string> utils::mesh_get_tokens(std::ifstream& mesh_file)
{
  std::string line;                // file line buffer
  std::vector<std::string> tokens; // tokenized line
  getline(mesh_file, line);
  boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
  return (tokens);
}

void utils::save_matrix(const std::string file_name, int bytes, char* data)
{
  std::ofstream data_file;
  data_file.open(file_name, std::ios::binary);
  data_file.write(data, bytes);
  data_file.close();
}
