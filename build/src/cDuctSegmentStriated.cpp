/*
 * cDuctSegmentStriated.cpp
 *
 *  Created on: 21/04/2021
 *      Author: jrugis
 */

//#include <thread>

#include "cCell.hpp"
#include "cDuctSegment.hpp"

cDuctSegmentStriated::cDuctSegmentStriated(cMiniGlandDuct* parent, int seg_number) : cDuctSegment(parent, seg_number){}

void cDuctSegmentStriated::step(){

  // concurrent cells step  
  //std::vector<std::thread> threads;
  //for(auto cell : cells) {
  //  threads.emplace_back([&](){cell->step();}); // NOTE: step function passed by reference
  //}
  //for(auto& t : threads) t.join(); // wait for all cell threads to complete

  //#pragma omp parallel num_threads(10)
  //{
  #pragma omp parallel for
  for(auto cell : cells){
	cell->step();
  }
   // for(unsigned int i=0; i<cells.size();i++) cells[i]->step();
  //}


  // combine cells fluid flow  --  TO DO
  // ....

  // NOTE - do this to the cells vector in the derived class access the the cell data 
  //  Base* basepointer = new Derived;
  //  static_cast<Derived*>(basepointer)->derived_int; // Can now access
  
  out << "<DuctSegmentStriated> step" << std::endl;
}
