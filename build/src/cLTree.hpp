/*
 * cLTree.hpp
 *
 *  Created on: 26/09/2021
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

  MatrixN3d nodes;     // lumen tree nodes
  MatrixN1d radii;     //            radi, per node
  MatrixN2i segs;      //            segments
  MatrixN1i types;     //            types, per segment ( 0 = duct, 1 = first acinus, 2 = second acinus, ...)

protected:

private:
  cMiniGland* parent;
  std::string id;
  std::ofstream out;     // runtime diagnostic file for this object
};

#endif /* CLTREE_H_ */
