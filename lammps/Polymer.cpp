// Polymer.cpp

#include <iostream>
#include <vector>
#include <memory>
#include <cstdlib>
#include <cmath>
#include <random>
#include <armadillo>
#include "Bead.hpp"
#include "Bond.hpp"
#include "Angle.hpp"
#include "Polymer.hpp"

using std::cout;
using std::endl;
using std::vector;
using std::shared_ptr;
using std::make_shared;
using namespace arma;

double getRand();
vec rosette(const double& r, const double& a, const double& k, 
	    const double& p, const double& t);
mat randRotation(double x1, double x2, double x3);

// Constructors
Polymer::Polymer() : Polymer {0, 1, 1, 1, false} {} 
Polymer::Polymer(int nBeads, int beadType, 
		 int bondType, int angleType, bool createBead){

  if (nBeads > 0){
    beads.reserve(nBeads);

    if (createBead){

      int i {};
      shared_ptr<Bead> bead {};

      bead = make_shared<Bead>();
      bead->setType(beadType);
      beads.push_back(bead);
      i++;
      
      if (i < nBeads){
	bead = make_shared<Bead>();
	bead->setType(beadType);
	beads.push_back(bead);
	beads[i-1]->addBondWith(bondType, beads[i]);
	i++;
      }
	
      for (; i < nBeads; i++){
	bead = make_shared<Bead>();
	bead->setType(beadType);
	beads.push_back(bead);
	beads[i-1]->addBondWith(bondType, beads[i]);
	beads[i-2]->addAngleWith(angleType, beads[i-1], beads[i]);
      }
    }
  }
}

// Accessor methods
shared_ptr<Bead> Polymer::getBead(int id){
  return beads[id];
}

vector<shared_ptr<Bead> >& Polymer::getBeads(){
  return beads;
}

int Polymer::getNumOfBeads(){
  return beads.size();
}

// Adding or removing beads
void Polymer::addBead(shared_ptr<Bead> bead){
  beads.push_back(bead);
}

void Polymer::addBead(int id, shared_ptr<Bead> bead){
  beads.insert(beads.begin()+id, bead);
}

void Polymer::addBead(int bondType, int angleType, shared_ptr<Bead> bead){
  addBead(getNumOfBeads(), bondType, angleType, bead);
}

void Polymer::addBead(int id, int bondType, int angleType,
					  shared_ptr<Bead> bead){

  int nbeads {getNumOfBeads()};

  // Check that the id is within the valid range
  if (id < 0 || id > nbeads)
	return;
  
  beads.insert(beads.begin()+id, bead);
  nbeads++;

  // Make sure the polymer remains connected by removing/adding bonds/angles
  if (id == 0){
	bead->addBondWith(bondType, beads[1]);
	bead->addAngleWith(angleType, beads[1], beads[2]);
  } else if (id == nbeads-1){
	bead->addBondWith(bondType, beads[nbeads-2]);
	bead->addAngleWith(angleType, beads[nbeads-2], beads[nbeads-3]);
  } else if (id == 1){
	bead->addBondWith(bondType, beads[0]);
	bead->addBondWith(bondType, beads[2]);
	bead->addAngleWith(angleType, beads[0], beads[2]);
	bead->addAngleWith(angleType, beads[2], beads[3]);
	beads[0]->removeBondWith(beads[2]);
	beads[0]->removeAngleWith(beads[2], beads[3]);
  } else if (id == nbeads-2){
	bead->addBondWith(bondType, beads[nbeads-1]);
	bead->addBondWith(bondType, beads[nbeads-3]);
	bead->addAngleWith(angleType, beads[nbeads-1], beads[nbeads-3]);
	bead->addAngleWith(angleType, beads[nbeads-3], beads[nbeads-4]);
	beads[nbeads-1]->removeBondWith(beads[nbeads-3]);
	beads[nbeads-1]->removeAngleWith(beads[nbeads-3], beads[nbeads-4]);
  } else if (id > 1 && id < nbeads-2){
	bead->addBondWith(bondType, beads[id-1]);
	bead->addBondWith(bondType, beads[id+1]);
	bead->addAngleWith(angleType, beads[id-2], beads[id-1]);
	bead->addAngleWith(angleType, beads[id-1], beads[id+1]);
	bead->addAngleWith(angleType, beads[id+1], beads[id+2]);
	beads[id-1]->removeBondWith(beads[id+1]);
	beads[id-2]->removeAngleWith(beads[id-1], beads[id+1]);
	beads[id-1]->removeAngleWith(beads[id+1], beads[id+2]);
  }
}

void Polymer::removeBead(int id){
  // Find the bond and angle type with neighbouring bead
  int bondType {1};
  int angleType {1};
  shared_ptr<Bond> bond {}; 
  shared_ptr<Angle> angle {};
  int bead1Index {id+1};
  int bead2Index {id+2};

  bond = beads[id]->getBondWith(beads[bead1Index]);
  angle = beads[id]->getAngleWith(beads[bead1Index], beads[bead2Index]);

  if (bond != nullptr)
	bondType = bond->getType();
  if (angle != nullptr)
	angleType = angle->getType();

  // Remove all bonds and angles
  beads[id]->removeAllBonds();
  beads[id]->removeAllAngles();
  
  // Make sure the polymer remains connected by adding bonds/angles
  beads[id-1]->addBondWith(bondType, beads[id+1]);
  beads[id-2]->addAngleWith(angleType, beads[id-1], beads[id+1]);
  beads[id-1]->addAngleWith(angleType, beads[id+1], beads[id+2]);

  // Erase the current bead
  beads.erase(beads.begin()+id);
}

void Polymer::removeAllBeads(){
  for (auto const& bead : beads){
    bead->removeAllBonds();
    bead->removeAllAngles();
  }
  beads.clear();
}

// Statistics of polymer
vector<double> Polymer::getCentreOfMass(double lx, double ly, double lz){
  double x {}, y {}, z {};
  for (auto const& b : beads){
    x += (b->getPosition(0) + lx*b->getBoundaryCount(0));
    y += (b->getPosition(1) + ly*b->getBoundaryCount(1));
    z += (b->getPosition(2) + lz*b->getBoundaryCount(2));
  }
  double numOfBeads = beads.size();
  x /= numOfBeads;
  y /= numOfBeads;
  z /= numOfBeads;
  return {x, y, z};
}

double Polymer::getGyrationRadius(double lx, double ly, double lz){
  vector<double> cm {getCentreOfMass(lx, ly, lz)};
  double dx {}, dy {}, dz {}, sum {};
  for (auto const& b : beads){
    dx = b->getPosition(0) + lx*b->getBoundaryCount(0) - cm[0];
    dy = b->getPosition(1) + ly*b->getBoundaryCount(1) - cm[1];
    dz = b->getPosition(2) + lz*b->getBoundaryCount(2) - cm[2];
    sum += dx*dx+dy*dy+dz*dz;
  }
  sum /= beads.size();
  sum = sqrt(sum);
  return sum;
}

// For handling bead, bond, angle listeners
void Polymer::addBeadListener(const shared_ptr<BeadListener>& l){
  for (auto const& b : beads){
    b->addBeadListener(l);
  }
}

void Polymer::removeBeadListener(const shared_ptr<BeadListener>& l){
  for (auto const& b : beads){
    b->removeBeadListener(l);
  }
}

void Polymer::addBondListener(const shared_ptr<BondListener>& l){
  for (auto const& b : beads){
    b->addBondListener(l);
  }
}

void Polymer::removeBondListener(const shared_ptr<BondListener>& l){
  for (auto const& b : beads){
    b->removeBondListener(l);
  }
}

void Polymer::addAngleListener(const shared_ptr<AngleListener>& l){
  for (auto const& b : beads){
    b->addAngleListener(l);
  }
}

void Polymer::removeAngleListener(const shared_ptr<AngleListener>& l){
  for (auto const& b : beads){
    b->removeAngleListener(l);
  }
}

// Static factory methods for creating polymers
shared_ptr<Polymer> 
Polymer::createRandomWalkPolymer(int nBeads, int beadType,
				 int bondType, int angleType,
				 double x0, double y0, double z0, 
				 double lx, double ly, double lz){
  // Initialise the random number generator
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> randDouble(0,1.0);
  double pi {M_PI};
  shared_ptr<Polymer> polymer {make_shared<Polymer>(nBeads, beadType,
						    bondType, angleType)};
  double x, y, z, r, costheta, sintheta, phi;
  shared_ptr<Bead> previous {};
  shared_ptr<Bead> current {};

  // Set the first bead to be centred at the origin
  previous = polymer->getBead(0);
  previous->setPosition(0, x0);
  previous->setPosition(1, y0);
  previous->setPosition(2, z0);

  for (int i {1}; i < nBeads; i++){
    current = polymer->getBead(i);
    do {
      r = randDouble(mt);
      costheta = 1.0-2.0*r;
      sintheta = sqrt(1-costheta*costheta);
      r = randDouble(mt);
      phi = 2.0*pi*r;
      x = previous->getPosition(0) + sintheta * cos(phi);
      y = previous->getPosition(1) + sintheta * sin(phi);
      z = previous->getPosition(2) + costheta;
    } while (fabs(x-x0) > lx/2.0 || fabs(y-y0) > ly/2.0 || fabs(z-z0) > lz/2.0);
    current->setPosition(0, x);
    current->setPosition(1, y);
    current->setPosition(2, z);
	previous = current;
  }
  return polymer;
}

shared_ptr<Polymer> 
Polymer::createRosettePolymer(int nBeads, int beadType, int bondType, 
			      int angleType, int beadsPerTurn, 
			      double r, double a, double k, double p,
			      double x0, double y0, double z0, 
			      double lx, double ly, double lz){
  // Initialise the random number generator
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> randDouble(0,1.0);
  double pi {M_PI};
  shared_ptr<Polymer> polymer {make_shared<Polymer>(nBeads, beadType,
						    bondType, angleType)};
  double tInc {beadsPerTurn / (2.0*pi)};
  
  // Determine height of rosette cylinder
  double height {nBeads / static_cast<double>(beadsPerTurn)};

  const vec centre {x0, y0, z0};

  // Generate random orientation
  bool outOfBound {true};
  double x1, x2, x3;
  mat rotate (3,3, fill::zeros);
  vec pos {0.0, 0.0, -height/2.0}; // Position in rotated frame
  double xlo {-lx/2.0}, xhi {lx/2.0};
  double ylo {-ly/2.0}, yhi {ly/2.0};
  double zlo {-lz/2.0}, zhi {lz/2.0};
  vec origin;
  
  while (outOfBound){
	outOfBound = false;
	// Gernerate random number
	x1 = randDouble(mt);
	x2 = randDouble(mt);
	x3 = randDouble(mt);

	// Generate rotation
	rotate = randRotation(x1, x2, x3);
	pos = rotate * pos; // Position in the lab frame
	origin = pos;
	
	// Check that the bottom rosette is within the box
	for (int i {}; i < beadsPerTurn; i++){
	  if (pos(0) < xlo || pos(0) > xhi ||
		  pos(1) < ylo || pos(1) > yhi ||
		  pos(2) < zlo || pos(2) > zhi){
		outOfBound = true;
		break;
	  }
	  pos = (rotate * rosette(r, a, k, p, i*tInc)) + origin;
	}
  }
  
  pos = origin;
  shared_ptr<Bead> bead {};

  for (int i {}; i < nBeads; i++){
	bead = polymer->getBead(i);
	// Set bead position
	bead->setPosition(0, pos(0));
	bead->setPosition(1, pos(1));
	bead->setPosition(2, pos(2));
	// Get next bead position
	pos = (rotate * rosette(r, a, k, p, i*tInc)) + origin;
  }
  return polymer;
}

double getRand(){
  return static_cast<double>(rand())/static_cast<double>(RAND_MAX);
}

// Define the rosette vector
vec rosette(const double& r, const double& a, const double& k, 
	    const double& p, const double& t){
  vec v (3);
  const double pi {M_PI};
  double cos2kt {cos(k*t)}; 
  cos2kt *= cos2kt;
  v(0) = r*(a+(1-a)*cos2kt*cos(t));
  v(1) = r*(a+(1-a)*cos2kt*sin(t));
  v(2) = p*t/2.0/pi;
  return v;
}


mat randRotation(double x1, double x2, double x3){
  const double pi2 {M_PI*2};
  double theta {x1 * pi2}; // Rotation about the pole (Z)
  double phi {x2 * pi2}; // Direction of pole deflection
  double z {x3 * 2.0}; // Magnitude of pole deflection

  double r {sqrt(z)};
  double vx {cos(phi) * r};
  double vy {sin(phi) * r};
  double vz {sqrt(2.0 - z)};

  double st {sin(theta)};
  double ct {cos(theta)};
  double sx {vx * ct + vy * st};
  double sy {vy * ct - vx * st};

  // Construct the rotation matrix
  mat rotation(3,3);
  rotation.at(0,0) = vx * sx - ct;
  rotation.at(0,1) = vx * sy + st;
  rotation.at(0,2) = vx * vz;
  
  rotation.at(1,0) = vy * sx - st;
  rotation.at(1,1) = vy * sy - ct;
  rotation.at(1,2) = vy * vz;
  
  rotation.at(2,0) = vz * sx;
  rotation.at(2,1) = vz * sy;
  rotation.at(2,2) = 1.0 - z;
  return rotation;
}
