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
//#include <algorithm>
//#include <execution>
#include <iomanip>
#include <iostream>
#include <string>
//#include <thread>

#include "global_defs.hpp"
#include "utils.hpp"

#include "cDuctSegment.hpp"
#include "cMiniGlandDuct.hpp"

cMiniGlandDuct::cMiniGlandDuct()
{
  id = "d1"; // there's only a single mini-gland duct
  out.open(id + DIAGNOSTIC_FILE_EXTENSION);
  utils::get_parameters(PARAMETER_FILE_NAME, p, out); // NOTE: all the parameters are in this file

  // get duct segment summary data
  std::string line;                    // file line buffer
  std::vector<std::string> tokens;     // tokenized line
  std::ifstream mesh_file(MESH_FILE_NAME); // open the mesh file
  if (not mesh_file.is_open()) { utils::fatal_error("mesh file " + std::string(MESH_FILE_NAME) + " could not be opened", out); }

  // get the duct node count
  while (getline(mesh_file, line)) {
  	boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
  	if (tokens[1]==std::string("duct_node") ) break;
  }
  int nnodes = atoi(tokens[2].c_str());

  // get the duct segment count
  while (getline(mesh_file, line)) {
  	boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
  	if (tokens[1]==std::string("duct_segment") ) break;
  }
  int nsegs = atoi(tokens[2].c_str());
  // skip over the rest of the header
  while (getline(mesh_file, line)) {
  	boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
  	if (tokens[0]==std::string("end_header") ) break;
  }
  // skip over the duct nodes
  for(int i=0; i<nnodes; i++) getline(mesh_file, line);	

  // get data for each duct segment and create duct segment objects
  for(int i=0; i<nsegs; i++){
	getline(mesh_file, line);
    boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
    int vertex_in = atoi(tokens[0].c_str());
    int vertex_out = atoi(tokens[1].c_str());
    int seg_type = atoi(tokens[4].c_str());
	if(seg_type == ACINUS) segments.push_back(new cDuctSegmentAcinus(this, i));
	else if(seg_type == INTERCALATED) segments.push_back(new cDuctSegmentIntercalated(this, i));
	else if(seg_type == STRIATED) segments.push_back(new cDuctSegmentStriated(this, i));
	seg_data.push_back(std::make_tuple(vertex_in, vertex_out, seg_type));
  }
    
  mesh_file.close();
  out << "<MiniGlandDuct> Segment count: " << segments.size() << std::endl;
}

cMiniGlandDuct::~cMiniGlandDuct() { 
  for(unsigned int i = 0; i<segments.size(); i++) delete segments[i]; // delete the duct segments
  out.close(); // close the diagnostic file
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

    // concurrent duct segment step  
    //std::vector<std::thread> threads;
    //for(auto seg : segments) {
    //  threads.emplace_back([&](){seg->step();}); // NOTE: step function passed by reference
    //}
    //for(auto& t : threads) t.join(); // wait for all duct segment threads to complete
    
	for(auto seg : segments) {
	  seg->step();
    }
	
	//std::for_each(std::execution::par_unseq, segments.begin(), segments.end(), [](auto&& seg)
	//{
	//  seg->step();
	//});


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
