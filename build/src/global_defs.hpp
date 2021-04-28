/*
 * global_defs.hpp
 *
 *  Created on: 21/04/2021
 *      Author: jrugis
 */

#ifndef DEFS_H_
#define DEFS_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>
//#include <cmath>
//#include <string>

#define PARAMETER_FILE_NAME "p.dat"
#define MESH_FILE_NAME "m.ply"
#define DIAGNOSTIC_FILE_EXTENSION ".out"

#define ACINUS 0
#define INTERCALATED 1
#define STRIATED 2

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
typedef Eigen::Vector3i Vector3i;
typedef Eigen::Vector4i Vector4i;

//************************************************************************
//************************************************************************
// cell mesh values
// struct sMeshVals {
//  int vertices_count;
//  int surface_triangles_count;
//  int tetrahedrons_count;
//  MatrixN3d vertices;          // 3x coordinate
//  MatrixN3i surface_triangles; // 3x vertex
//  MatrixN4i tetrahedrons;      // 4x vertex
//};

//************************************************************************
// thermodynamic constants
//#define R 8.314462100000000
//#define T 310
//#define F 9.645833650000000e4
// const double RTF = 1000.0 * R * T / F;

//************************************************************************
//************************************************************************

#endif /* DEFS_H_ */
