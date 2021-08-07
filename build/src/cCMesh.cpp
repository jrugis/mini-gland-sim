/*
 * cCMesh.cpp
 *
 *  Created on: 26/07/2021
 *      Author: jrugis
 */

#include <iostream>
#include <string>

#include "global_defs.hpp"
#include "utils.hpp"

#include "cCMesh.hpp"

cCMesh::cCMesh(int cell_number, std::ofstream& out)
{
  // get cell mesh data
  std::ifstream mesh_file;
  utils::mesh_open(mesh_file, out);                                             // open the mesh file
  int nacinii = utils::mesh_get_count(mesh_file, std::string("acinii"));        // number of acinii
  int nducts = utils::mesh_get_count(mesh_file, std::string("duct"));           //           ducts
  int nlnodes = utils::mesh_get_count(mesh_file, std::string("lumen_node"));    //         lumen nodes
  int nlsegs = utils::mesh_get_count(mesh_file, std::string("lumen_segment"));  //       lumen segments
  int tncells = utils::mesh_get_count(mesh_file, std::string("cell"));       // get the total mesh cell count
  int tnverts = utils::mesh_get_count(mesh_file, std::string("vertex"));     // get the total mesh vertex count
  int tnfaces = utils::mesh_get_count(mesh_file, std::string("face"));       // get the total mesh face count
  utils::mesh_end_header(mesh_file);                                         // skip over the rest of the header

  utils::mesh_skip_lines(mesh_file, nacinii + nducts + nlnodes + nlsegs + cell_number);     // skip to the correct cell
  std::vector<std::string> tokens = utils::mesh_get_tokens(mesh_file); // get cell data
  nverts = std::stoi(tokens[0]);                                       // number of vertices
  int iverts = std::stoi(tokens[1]);                                   //    vertices start index
  nfaces = std::stoi(tokens[2]);                                       // number of faces
  int ifaces = std::stoi(tokens[3]);                                   //    faces start index
  ntets = std::stoi(tokens[4]);                                        // mumber of tetrahedrons
  int itets = std::stoi(tokens[5]);                                    //    tetrahedrons start index
  utils::mesh_skip_lines(mesh_file, tncells - cell_number - 1);        // skip over the rest of cells

  out << "<CMesh> cell_number: " << cell_number+1 << std::endl;
  out << "<CMesh> number of vertices: " << nverts << std::endl;
  out << "<CMesh> number of faces: " << nfaces << std::endl;
  out << "<CMesh> number of tetrahedrons: " << ntets << std::endl;

  // get vertices for this cell
  utils::mesh_skip_lines(mesh_file, iverts); // skip to the correct vertex block
  verts.resize(nverts, Eigen::NoChange);     // resize the verts matrix
  for (int i = 0; i < nverts; i++) {         // get vertices
    std::vector<std::string> tokens = utils::mesh_get_tokens(mesh_file);
	verts.row(i) = Vector3d(std::stod(tokens[0]), std::stod(tokens[1]), std::stod(tokens[2])); // save vertex
  }
  utils::mesh_skip_lines(mesh_file, tnverts - iverts - nverts); // skip over the remaining verts

  // get face data for this cell
  utils::mesh_skip_lines(mesh_file, ifaces);  // skip to the correct faces block
  faces.resize(nfaces, Eigen::NoChange);      // resize the faces matrix
  face_types.resize(nfaces, Eigen::NoChange); // resize the face types matrix
  for (int i = 0; i < nfaces; i++) {          // get vertex indices and face type
	std::vector<std::string> tokens = utils::mesh_get_tokens(mesh_file);
	faces.row(i) = Vector3i(std::stoi(tokens[0]) - iverts, std::stoi(tokens[1]) - iverts,
	               std::stoi(tokens[2]) - iverts); // save vertex indices
  }
  utils::mesh_skip_lines(mesh_file, tnfaces - ifaces - nfaces); // skip over the remaining faces
  
  // get tetrahedron data for this cell
  if (ntets > 0) {                                                         // NOTE: only acinus cells have tets!
	utils::mesh_skip_lines(mesh_file, itets);                              // skip to the correct tets block
    tets.resize(ntets, Eigen::NoChange);                                   // resize the tets matrix
	for (int i = 0; i < ntets; i++) {                                      // get vertex indices
	  std::vector<std::string> tokens = utils::mesh_get_tokens(mesh_file); // four vertex indices
	  tets.row(i) = Vector4i(std::stoi(tokens[0]) - itets, std::stoi(tokens[1]) - itets, std::stoi(tokens[2]) - itets,
	                std::stoi(tokens[3]) - itets); // save them
	}
  }
  mesh_file.close();

  // calculate all of the face centers
  utils::calc_tri_centers(face_centers, verts, faces);
  
  // calculate all of the face areas
  utils::calc_tri_areas(face_areas, verts, faces);
}

cCMesh::~cCMesh()
{
}
