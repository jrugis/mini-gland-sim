/*
 * cMiniGlandDuct.cpp
 *
 *  Created on: 21/04/2021
 *      Author: jrugis
 */

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <string>

#include "global_defs.hpp"
#include "utils.hpp"

#include "cDuctSegment.hpp"
#include "cMiniGlandDuct.hpp"

cMiniGlandDuct::cMiniGlandDuct()
{
  id = "d1"; // there's only a single mini-gland duct
  out.open(id + DIAGNOSTIC_FILE_EXTENSION);
  utils::get_parameters(PARAMETER_FILE_NAME, p, out); // NOTE: all the parameters are in this file

  // TO DO - get segment data and create segments
  // push_back the different seg types
  
  // create duct segments
  //int nsegs = get_element_count();
  //for(int i =0; i<nsegs; i++){
  //	  segments.push_back(new cDuctSegment(this, i));
  //}
  //out << "<MiniGlandDuct> Segment count: " << segments.size() << std::endl;
}

cMiniGlandDuct::~cMiniGlandDuct() { 
  for(int i = 0; i<segments.size(); i++) delete segments[i]; // delete the duct segments
  out.close(); // close the diagnostic file
}

std::tuple<int, std::vector<int>> cMiniGlandDuct::get_segment_data()
{
  std::string line;                    // file line buffer
  std::vector<std::string> tokens;     // tokenized line
  std::ifstream mesh_file(MESH_FILE_NAME); // open the mesh file
  if (not mesh_file.is_open()) { utils::fatal_error("cell file " + std::string(MESH_FILE_NAME) + " could not be opened", out); }
  while (getline(mesh_file, line)) {
	boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
	if (tokens[1]=="duct_segment") break;
  }
  mesh_file.close();
  return(atoi(tokens[2].c_str()));
}

void cMiniGlandDuct::run()
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
