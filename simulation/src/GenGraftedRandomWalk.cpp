/* GenGraftedRandomWalk.cpp
 * This is a code that generaetes a random walk polymer.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <memory>
#include <cstdlib>
#include <cmath>
#include "Bead.hpp"
#include "Polymer.hpp"
#include "LAMMPS.hpp"
#include "DataManager.hpp"

using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::string;
using std::shared_ptr;
using std::make_shared;

int getBeadType(double fracOfLam, double fracOfHet);

int main(int argc, char * argv[]){
  if (argc != 8){
    cout << "Usage: [num of beads] [lx] [ly] [lz] [buffer] " 
	 << "[output file] [output map file]" << endl;
    return 1;
  }

  int argi {};
  int numOfBeads {stoi(string(argv[++argi]), nullptr, 10)};
  double lx {stod(string(argv[++argi]), nullptr)};
  double ly {stod(string(argv[++argi]), nullptr)};
  double lz {stod(string(argv[++argi]), nullptr)};
  double buffer {stod(string(argv[++argi]), nullptr)};
  string outFile (argv[++argi]);
  string outMapFile (argv[++argi]);

  // Generate polymer
  shared_ptr<LAMMPS> lammps = make_shared<LAMMPS>(lx, ly, lz);
  shared_ptr<Polymer> polymer {};
  shared_ptr<Bead> bead, head {};
  
  lammps->setTypesOfBeads(2);
  lammps->setTypesOfBonds(2);
  lammps->setTypesOfAngles(1);

  polymer = lammps->createRandomWalkPolymer(0, numOfBeads, 1, 1, 1, 
					    0.0, 0.0, 0.0, 
					    lx-buffer, ly-buffer, lz-buffer);
  for (int i {}; i < numOfBeads; i++){
    bead = lammps->getPolymer(0)->getBead(i);
    bead->setLabel(1);
    bead->setType(1);
  }

  // Add an extra bead to tether the beginning end of the polymer
  // to the adsorptive wall
  // Add the tethering bead at the centre of the wall
  bead = lammps->createBead(1, 0.0, 0.0, lz/2.0-0.01, 0.0, 0.0, 0.0, 0, 0, 0, 2, 1);
  head = lammps->getPolymer(0)->getBead(0);
  bead->addBondWith(2, head);

  // Write the input file
  lammps->exportData(outFile, outMapFile);
}
