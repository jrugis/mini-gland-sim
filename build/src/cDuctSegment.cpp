/*
 * cDuctSegment.cpp
 *
 *  Created on: 21/04/2021
 *      Author: jrugis
 */

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/tokenizer.hpp>

//#include <iomanip>
#include <iostream>
#include <string>

#include "global_defs.hpp"
#include "utils.hpp"

#include "cCell.hpp"
#include "cMiniGlandDuct.hpp"
#include "cDuctSegment.hpp"

cDuctSegment::cDuctSegment(cMiniGlandDuct* _parent, int _seg_number) : parent(_parent), seg_number(_seg_number)
{

  id = parent->id + "s" + std::to_string(seg_number+1); // NOTE: one based segment id
  out.open(id + DIAGNOSTIC_FILE_EXTENSION);

  utils::get_parameters(PARAMETER_FILE_NAME, p, out); // NOTE: all the parameters are in this file
  get_segment_data();
  out << "<DuctSegment> Cell count: " << cells.size() << std::endl;
}

cDuctSegment::~cDuctSegment() { 
  for(int i = 0; i<cells.size(); i++) delete cells[i]; // delete the cells
  out.close(); 
}

void cDuctSegment::get_segment_data(){
  std::string line;                    // file line buffer
  std::vector<std::string> tokens;     // tokenized line
  std::vector<Vector3d> verts;

  std::ifstream mesh_file(MESH_FILE_NAME); // open the mesh file
  if (not mesh_file.is_open()) { utils::fatal_error("cell file " + std::string(MESH_FILE_NAME) + " could not be opened", out); }

  // get the duct node vertex count
  while (getline(mesh_file, line)) {
  	boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
  	if (tokens[1]=="duct_node" ) break;
  }
  int nnodes = atoi(tokens[2].c_str());

  // skip over the rest of the header
  while (getline(mesh_file, line)) {
  	boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
  	if (tokens[0]=="end_header" ) break;
  }

  // get the node data
  for(int i=0; i<nnodes; i++){
  	getline(mesh_file, line);
	boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
	verts.push_back(Vector3d(atof(tokens[0].c_str()), atof(tokens[1].c_str()), atof(tokens[2].c_str())));
  }

  // skip to the correct segment and get the data
  for(int i=0; i<seg_number+1; i++) getline(mesh_file, line);
  boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
  vertex_in = verts[atoi(tokens[0].c_str())];
  vertex_out = verts[atoi(tokens[1].c_str())];
  inner_diameter = atof(tokens[2].c_str());
  outer_diameter = atof(tokens[3].c_str());
  seg_type = atoi(tokens[4].c_str());

ToDO check this...  
  // create the cells for this segment
  int ncells = atoi(tokens[5].c_str());
  for(int i =0; i<ncells; i++){
	 cells.push_back(new cCell(this, atoi(tokens[6+i].c_str())));
  }

  mesh_file.close();
}

void cDuctSegment::run()
{
  /*  double t = 0.0;
    double solver_dt = p.at("delT");
    double error;
    struct timespec start, end;
    double elapsed;

    // simulation time stepping and synchronization
    clock_gettime(CLOCK_REALTIME, &start);
    while ((p.at("totalT") - t) > 0.000001) { // HARD CODED: assumes solver_dt always > 1us
      error = snd_recv(t, solver_dt);         // invoke the calcium solver
      if (error != 0.0) {
        // ...   change time step???
      }
      // invoke lumen step ???

      clock_gettime(CLOCK_REALTIME, &end);
      elapsed = (end.tv_sec - start.tv_sec) + ((end.tv_nsec - start.tv_nsec) / 1000000000.0);
      out << std::fixed << std::setprecision(3);
      out << "<Acinus> step duration: " << elapsed << "s" << std::endl;
      start = end;
      t += solver_dt;
    }

    // instruct cells to finish
    snd_recv(t, 0.0);
    */
}
