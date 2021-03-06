/*
 * global_defs.hpp
 *
 *  Created on: 21/04/2021
 *      Author: jrugis
 */

#ifndef DEFS_H_
#define DEFS_H_

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
//#include <cmath>
#include <string>
#include <limits>

#define PARAMETER_FILE_NAME "p.ini"
#define DUCT_FILE_NAME "./meshes/_duct.ply"
#define MESH_FILE_DIR "./meshes"
#define DIAGNOSTIC_FILE_EXTENSION ".out"

enum face_types{APICAL, BASOLATERAL, BASAL};                  // face types
enum cell_groups{STRIATED, INTERCALATED, ACINAR_1, ACINAR_2}; // cell groups
enum lt_groups{DUCT, ACINUS_1, ACINUS_2};       // lumen tree segment groups 

//************************************************************************
// some convenience typedefs
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> MatrixN1d;
typedef Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajorBit> MatrixN3d;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixNNd;
typedef Eigen::SparseMatrix<double> SparceMatrixd;

typedef Eigen::Matrix<int, Eigen::Dynamic, 1> MatrixN1i;
typedef Eigen::Matrix<int, Eigen::Dynamic, 2, Eigen::RowMajorBit> MatrixN2i;
typedef Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajorBit> MatrixN3i;
typedef Eigen::Matrix<int, Eigen::Dynamic, 4, Eigen::RowMajorBit> MatrixN4i;
typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> MatrixNNi;

typedef Eigen::Vector3d Vector3d;
typedef Eigen::Vector2i Vector2i;
typedef Eigen::Vector3i Vector3i;
typedef Eigen::Vector4i Vector4i;

typedef Eigen::Array<double, 1, Eigen::Dynamic> Array1Nd;

//************************************************************************
// thermodynamic constants
#define R 8.13144621
#define T 310
#define F_const 96485.3329
#define V_w 18e12
const double RTF = 1000.0 * R * T / F_const;

//************************************************************************
// parameter constants
const double PARAMETER_DEFAULT_REAL = std::numeric_limits<double>::lowest();
const int PARAMETER_DEFAULT_INTEGER = std::numeric_limits<int>::lowest();
const std::string PARAMETER_DEFAULT_STRING = "<UNSET>";

//************************************************************************
//************************************************************************

#endif /* DEFS_H_ */
