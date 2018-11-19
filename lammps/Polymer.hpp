// Polymer.hpp

#ifndef POLYMER_HPP
#define POLYMER_HPP

#include <vector>
#include <memory>
#include "Bead.hpp"

using std::vector;
using std::shared_ptr;

class Polymer {

private:
  vector<shared_ptr<Bead> > beads {};

public:
  // Constructor
  Polymer(); 
  Polymer(int numOfBeads, int beadType, 
	  int bondType, int angleType, bool createBead = true);

  // Accessor methods
  shared_ptr<Bead> getBead(int id);
  vector<shared_ptr<Bead> >& getBeads();

  int getNumOfBeads();

  // For changing the polymer
  void addBead(shared_ptr<Bead> bead);
  void addBead(int id, shared_ptr<Bead> bead);
  void addBead(int bondType, int angleType, shared_ptr<Bead> bead);
  void addBead(int id, int bondType, int angleType, shared_ptr<Bead> bead);
  void removeBead(int id);
  void removeAllBeads();
  
  // Statistics of polymer
  vector<double> getCentreOfMass(double lx, double ly, double lz);
  double getGyrationRadius(double lx, double ly, double lz);

  // For handling bead, bond, angle listeners
  void addBeadListener(const shared_ptr<BeadListener>& listener);
  void removeBeadListener(const shared_ptr<BeadListener>& listener);
  void addBondListener(const shared_ptr<BondListener>& listener);
  void removeBondListener(const shared_ptr<BondListener>& listener);
  void addAngleListener(const shared_ptr<AngleListener>& listener);
  void removeAngleListener(const shared_ptr<AngleListener>& listener);
  
  // Static factory methods
  static shared_ptr<Polymer> 
  createRandomWalkPolymer(int numOfBeads, int beadType, 
			  int bondType, int angleType,
			  double x, double y, double z,
			  double lx, double ly, double lz);

  static shared_ptr<Polymer> 
  createRosettePolymer(int nBeads, int beadType, 
		       int bondType, int angleType, int beadsPerTurn, 
		       double r, double a, double k, double p,
		       double x0, double y0, double z0, 
		       double lx, double ly, double lz);
  };

#endif
