/*
 * cSIMesh.cpp
 *
 *  Created on: 26/09/2021
 *      Author: jrugis
 */

#include <iostream>
#include <string>

#include "global_defs.hpp"
#include "utils.hpp"

#include "cSIMesh.hpp"

cSIMesh::cSIMesh(std::string fname, std::ofstream& out)
{
  std::ifstream mesh_file;
  utils::mesh_open(mesh_file, fname, out);                 // open the mesh file
  int nverts = utils::mesh_get_count(mesh_file, "vertex"); // get mesh vertex count
  int nfaces = utils::mesh_get_count(mesh_file, "face");   // get mesh face count
  utils::mesh_end_header(mesh_file);                       // skip over the rest of the header
  out << "<CMesh> number of vertices: " << nverts << std::endl;
  out << "<CMesh> number of faces: " << nfaces << std::endl;

  verts.resize(nverts, Eigen::NoChange);
  for (int i = 0; i < nverts; i++) {         // get vertices
	std::vector<std::string> tokens = utils::mesh_get_tokens(mesh_file);
	verts.row(i) = Vector3d(std::stod(tokens[0]), std::stod(tokens[1]), std::stod(tokens[2]));
  }
  faces.resize(nfaces, Eigen::NoChange);
  face_types.resize(nfaces, Eigen::NoChange);
  duct_idx.resize(nfaces, Eigen::NoChange);
  dist_from_duct.resize(nfaces, Eigen::NoChange);
  dist_along_duct.resize(nfaces, Eigen::NoChange);
  for (int i = 0; i < nfaces; i++) {         // get faces and face related data
	std::vector<std::string> tokens = utils::mesh_get_tokens(mesh_file);
	faces.row(i) = Vector3i(std::stoi(tokens[1]), std::stoi(tokens[2]), std::stoi(tokens[3]));
    face_types(i) = std::stoi(tokens[4]);
    duct_idx(i) = std::stoi(tokens[5]);
    dist_from_duct(i) = std::stod(tokens[6]);
    dist_along_duct(i) = std::stod(tokens[7]);
  }
  mesh_file.close();
  utils::calc_tri_areas(face_areas, verts, faces); // calculate all of the face areas
}

cSIMesh::~cSIMesh()
{
}
