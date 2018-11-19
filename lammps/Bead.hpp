// Bead.h

#ifndef BEAD_HPP
#define BEAD_HPP

#include <set>
#include <vector>
#include <memory>
#include "Bond.hpp"
#include "Angle.hpp"
#include "BeadListener.hpp"
#include "BondListener.hpp"
#include "AngleListener.hpp"

using std::set;
using std::vector;
using std::weak_ptr;
using std::shared_ptr;

class Bead : public std::enable_shared_from_this<Bead> {
  
private:
  double position [3];
  double velocity [3];
  int boundaryCount [3];
  int type;
  int label;

  vector<shared_ptr<Bond> > bondList {};
  vector<shared_ptr<Angle> > angleList {};
  
  vector<weak_ptr<BeadListener> > beadListeners {};
  vector<weak_ptr<BondListener> > bondListeners {};
  vector<weak_ptr<AngleListener> > angleListeners {};

  void addBond(shared_ptr<Bond> bond);
  void addAngle(shared_ptr<Angle> angle);

  void notifyBeadTypeChange(int oldType, int newType);
  void notifyBeadLabelChange(int oldLabel, int newLabel);
  void notifyBondCreation(const shared_ptr<Bond>& bond);
  void notifyBondRemoval(const shared_ptr<Bond>& bond);
  void notifyAngleCreation(const shared_ptr<Angle>& angle);
  void notifyAngleRemoval(const shared_ptr<Angle>& angle);

  
public:
  // Constructors
  Bead(double x, double y, double z,
       double vx, double vy, double vz,
       int nx, int ny, int nz, int type, int label);
  Bead(double x, double y, double z);
  Bead();
  
  // Accesor methods
  void setPosition(int dim, double value);
  double getPosition(int dim);
  void setVelocity(int dim, double value);
  double getVelocity(int dim);
  void setBoundaryCount(int dim, int count);
  int getBoundaryCount(int dim);
  void setType(int type);
  int getType();
  void setLabel(int label);
  int getLabel();
  int getNumOfBonds();
  int getNumOfAngles();
  shared_ptr<Bond> getBondWith(const shared_ptr<Bead>& bead);
  shared_ptr<Angle> getAngleWith(const shared_ptr<Bead>& bead1, 
				 const shared_ptr<Bead>& bead2);
  vector<shared_ptr<Bond> >& getBonds();
  vector<shared_ptr<Angle> >& getAngles();

  // For modifying bonds and angles
  void addBondWith(int type, const shared_ptr<Bead>& bead);
  void removeBond(const shared_ptr<Bond>& bond);
  void removeBondWith(const shared_ptr<Bead>& bead);
  void addAngleWith(int type, const shared_ptr<Bead>& bead1, 
		    const shared_ptr<Bead>& bead2);
  void removeAngle(const shared_ptr<Angle>& angle);
  void removeAngleWith(const shared_ptr<Bead>& bead1, 
		       const shared_ptr<Bead>& bead2);
  void removeAllBonds();
  void removeAllAngles();

  // For handling bead listeners
  void addBeadListener(const shared_ptr<BeadListener>& listener);
  void removeBeadListener(const shared_ptr<BeadListener>& listener);
  void addBondListener(const shared_ptr<BondListener>& listener);
  void removeBondListener(const shared_ptr<BondListener>& listener);
  void addAngleListener(const shared_ptr<AngleListener>& listener);
  void removeAngleListener(const shared_ptr<AngleListener>& listener);

};

#endif
