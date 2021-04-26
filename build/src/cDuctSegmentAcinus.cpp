/*
 * cDuctSegmentAcinus.cpp
 *
 *  Created on: 21/04/2021
 *      Author: jrugis
 */

#include "cDuctSegment.hpp"

cDuctSegmentAcinus::cDuctSegmentAcinus(cMiniGlandDuct* parent, int seg_number) : cDuctSegment(parent, seg_number){}

void cDuctSegmentAcinus::step(){
  // invoke cells step
  // for(auto cell : cells){
  //   cell->cell_step();		
  // }
  // combine cells fluid flow  --  TO DO
  // NOTE - do this to the cells vector in the derived class access the the cell data 
  //  Base* basepointer = new Derived;
  //  static_cast<Derived*>(basepointer)->derived_int; // Can now access
  // ....
  
  out << "<DuctSegmentAcinus> step" << std::endl;

}