/*
 * cMiniGland.cpp
 *
 *  Created on: 26/09/2021
 *      Author: jrugis
 */

#include <iomanip>
#include <iostream>
#include <string>

#include "global_defs.hpp"
#include "utils.hpp"

#include "cAcinus.hpp"
#include "cDuct.hpp"
#include "cLTree.hpp"
#include "cMiniGland.hpp"

cMiniGland::cMiniGland()
{
  id = "_mini-gland"; // there's only one mini-gland object
  out.open(id + DIAGNOSTIC_FILE_EXTENSION);
  p = new INIReader(PARAMETER_FILE_NAME);
  if (p->ParseError() < 0) {
    utils::fatal_error("Error opening parameter file", out);
  }
  ltree = new cLTree(this);  // create lumen tree object
  duct = new cDuct(this);    // create duct object
}

cMiniGland::~cMiniGland()
{
  delete ltree;
  delete duct;
  //for (unsigned int i = 0; i < acinii.size(); i++) delete acinii[i]; // delete the acinii
  out.close();     // close the diagnostic file
  delete p;
}

void cMiniGland::run()
{
  double t = 0.0;
  double solver_dt = utils::get_parameter_real(p, "time", "delT", out);
  double totalT = utils::get_parameter_real(p, "time", "totalT", out);
  struct timespec start, end;
  double elapsed;

  // simulation time stepping and synchronization
  clock_gettime(CLOCK_REALTIME, &start);
  while ((totalT - t) > 0.000001) { // HARD CODED: assumes solver_dt always > 1us

	// in parallel?
    // solve acinii fluid flow...
    
    // feed fluid flow into duct and 
	//    solve duct fluid flow...
	duct->step(t, solver_dt);

    // output step running time
    clock_gettime(CLOCK_REALTIME, &end);
    elapsed = (end.tv_sec - start.tv_sec) + ((end.tv_nsec - start.tv_nsec) / 1000000000.0);
    out << std::fixed << std::setprecision(3);
    out << "<MiniGland> current time: " << t + solver_dt << "; step duration: " << elapsed << "s" << std::endl;
    start = end;
    t += solver_dt;
  }
}
