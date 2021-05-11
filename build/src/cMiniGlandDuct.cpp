/*
 * cMiniGlandDuct.cpp
 *
 *  Created on: 21/04/2021
 *      Author: jrugis
 */

#include <iomanip>
#include <iostream>
#include <string>

#include "global_defs.hpp"
#include "utils.hpp"

#include "cDuctSegment.hpp"
#include "cDuctSegmentStriated.hpp"
#include "cMiniGlandDuct.hpp"

cMiniGlandDuct::cMiniGlandDuct()
{
  id = "d1"; // there's only a single mini-gland duct
  out.open(id + DIAGNOSTIC_FILE_EXTENSION);
  utils::get_parameters(PARAMETER_FILE_NAME, p, out); // NOTE: all the parameters are in this file
  
  // get duct segment mesh data
  std::ifstream mesh_file;
  utils::mesh_open(mesh_file, out);                                          // open the mesh file
  int nnodes = utils::mesh_get_count(mesh_file, std::string("duct_node"));   // get the duct node count
  int nsegs = utils::mesh_get_count(mesh_file, std::string("duct_segment")); // get the duct segment count
  utils::mesh_end_header(mesh_file);                                         // skip over the rest of the header
  utils::mesh_skip_lines(mesh_file, nnodes);                                 // skip over the duct node data
  
  // get data for each duct segment and create duct segment objects
  for (int i = 0; i < nsegs; i++) {
    std::vector<std::string> tokens = utils::mesh_get_tokens(mesh_file);
    int vertex_in = std::stoi(tokens[0]);
    int vertex_out = std::stoi(tokens[1]);
    int seg_type = std::stoi(tokens[4]);
    if (seg_type == ACINUS)
	  segments.push_back(new cDuctSegmentAcinus(this, i));
    else if (seg_type == INTERCALATED)
	  segments.push_back(new cDuctSegmentIntercalated(this, i));
    else if (seg_type == STRIATED)
	  segments.push_back(new cDuctSegmentStriated(this, i));
    seg_data.push_back(std::make_tuple(vertex_in, vertex_out, seg_type));
  }

  mesh_file.close();
  out << "<MiniGlandDuct> Segment count: " << segments.size() << std::endl;
}

cMiniGlandDuct::~cMiniGlandDuct()
{
  for (unsigned int i = 0; i < segments.size(); i++) delete segments[i]; // delete the duct segments
  out.close();                                                           // close the diagnostic file
}

void cMiniGlandDuct::run()
{
  double t = 0.0;
  double solver_dt = p.at("delT");
  struct timespec start, end;
  double elapsed;

  // simulation time stepping and synchronization
  clock_gettime(CLOCK_REALTIME, &start);
  while ((p.at("totalT") - t) > 0.000001) { // HARD CODED: assumes solver_dt always > 1us

#pragma omp parallel for
	for (auto seg : segments) { seg->step(); }

    // combine duct segment fluid flow  --  TO DO
    // ....

    // output step running time
    clock_gettime(CLOCK_REALTIME, &end);
    elapsed = (end.tv_sec - start.tv_sec) + ((end.tv_nsec - start.tv_nsec) / 1000000000.0);
    out << std::fixed << std::setprecision(3);
    out << "<MiniGlandDuct> step duration: " << elapsed << "s" << std::endl;
    start = end;
    t += solver_dt;
  }
}
