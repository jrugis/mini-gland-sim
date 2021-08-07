/*
 * cLTree.hpp
 *
 *  Created on: 26/07/2021
 *      Author: jrugis
 */

#ifndef CLTREE_H_
#define CLTREE_H_

#include <fstream>
#include <string>
#include <vector>

#include "global_defs.hpp"

class cMiniGland;

class cLTree {
public:
  cLTree(cMiniGland* parent);
  ~cLTree();

  MatrixN3d nodes;     // the lumen tree nodes
  MatrixN1d radii;     //                diameters
  MatrixN2i segs;      //                segments

protected:

private:
  cMiniGland* parent;
  std::string id;
  std::ofstream out;     // runtime diagnostic file for this object
};

#endif /* CLTREE_H_ */
