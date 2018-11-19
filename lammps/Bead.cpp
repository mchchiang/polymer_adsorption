// Bead.cpp

#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>
#include "Bead.hpp"
#include "Bond.hpp"
#include "Angle.hpp"
#include "BeadListener.hpp"
#include "BondListener.hpp"
#include "AngleListener.hpp"

using std::cout;
using std::endl;
using std::pair;
using std::weak_ptr;
using std::shared_ptr;
using std::make_shared;

// Constructors
Bead::Bead(double x, double y, double z,
		   double vx, double vy, double vz,
		   int nx, int ny, int nz, int t, int l) :
  position {x, y, z}, velocity {vx, vy, vz},
  boundaryCount {nx, ny, nz}, type {t}, label {l} {
    bondList.reserve(0);
    angleList.reserve(0);
}

Bead::Bead(double x, double y, double z) :
  Bead {x, y, z, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0} {}

Bead::Bead() :
  Bead {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0} {}

// Accessor methods
void Bead::setPosition(int dim, double value){
  position[dim] = value;
}

double Bead::getPosition(int dim){
  return position[dim];
}

void Bead::setVelocity(int dim, double value){
  velocity[dim] = value;
}

double Bead::getVelocity(int dim){
  return velocity[dim];
}

void Bead::setBoundaryCount(int dim, int count){
  boundaryCount[dim] = count;
}

int Bead::getBoundaryCount(int dim){
  return boundaryCount[dim];
}

void Bead::setType(int t){
  int oldType {type};
  type = t;
  if (oldType != type){
    notifyBeadTypeChange(oldType, type);
  }
}

int Bead::getType(){
  return type;
}

void Bead::setLabel(int l){
  int oldLabel {label};
  label = l;
  if (oldLabel != label){
    notifyBeadLabelChange(oldLabel, label);
  }
}

int Bead::getLabel(){
  return label;
}

int Bead::getNumOfBonds(){
  return bondList.size();
}

int Bead::getNumOfAngles(){
  return angleList.size();
}

shared_ptr<Bond> Bead::getBondWith(const shared_ptr<Bead>& bead){
  for (auto it = bondList.begin(); it != bondList.end(); it++){
    shared_ptr<Bond> bond = *it;
    if (bond->getBead(0) == bead || bond->getBead(1) == bead){
      return bond;
    }
  }
  return nullptr;
}

shared_ptr<Angle> Bead::getAngleWith(const shared_ptr<Bead>& bead1,
				     const shared_ptr<Bead>& bead2){
  for (auto it = angleList.begin(); it != angleList.end(); it++){
    shared_ptr<Angle> angle = *it;
    shared_ptr<Bead> b1 = angle->getBead(0);
    shared_ptr<Bead> b2 = angle->getBead(1);
    shared_ptr<Bead> b3 = angle->getBead(2);
    if ((bead1 == b1 && (bead2 == b2 || bead2 == b3)) ||
	(bead1 == b2 && (bead2 == b1 || bead2 == b3)) ||
	(bead1 == b3 && (bead2 == b1 || bead2 == b2))){
      return angle;
    }
  }
  return nullptr;
}

vector<shared_ptr<Bond> >& Bead::getBonds(){
  return bondList;
}

vector<shared_ptr<Angle> >& Bead::getAngles(){
  return angleList;
}

void Bead::addBond(shared_ptr<Bond> bond){
  bondList.push_back(bond);
}

void Bead::addBondWith(int t, const shared_ptr<Bead>& bead){
  shared_ptr<Bond> bond {make_shared<Bond>(t, shared_from_this(),bead)};
  bead->addBond(bond);
  addBond(bond);
  notifyBondCreation(bond);
}

void Bead::removeBond(const shared_ptr<Bond>& bond){
  auto it = std::find(bondList.begin(), bondList.end(), bond);
  if (it != bondList.end()){
    bondList.erase(it);
  }
}

void Bead::removeBondWith(const shared_ptr<Bead>& bead){
  for (auto it = bondList.begin(); it != bondList.end(); it++){
    shared_ptr<Bond> bond = *it;
    if (bond->getBead(0) == bead || bond->getBead(1) == bead){
      notifyBondRemoval(bond);
      bead->removeBond(bond);
      bondList.erase(it);
      break;
    }
  }
}

void Bead::addAngle(shared_ptr<Angle> angle){
  angleList.push_back(angle);
}

void Bead::addAngleWith(int t, const shared_ptr<Bead>& bead1, 
			const shared_ptr<Bead>& bead2){
  shared_ptr<Angle> angle {
    make_shared<Angle>(t, shared_from_this(), bead1, bead2)};
  bead1->addAngle(angle);
  bead2->addAngle(angle);
  addAngle(angle);
  notifyAngleCreation(angle);
}

void Bead::removeAngle(const shared_ptr<Angle>& angle){
  auto it = std::find(angleList.begin(), angleList.end(), angle);
  if (it != angleList.end()){
    angleList.erase(it);
  }
}

void Bead::removeAngleWith(const shared_ptr<Bead>& bead1, 
			   const shared_ptr<Bead>& bead2){
  for (auto it = angleList.begin(); it != angleList.end(); it++){
    shared_ptr<Angle> angle = *it;
    shared_ptr<Bead> b1 = angle->getBead(0);
    shared_ptr<Bead> b2 = angle->getBead(1);
    shared_ptr<Bead> b3 = angle->getBead(2);
    if ((bead1 == b1 && (bead2 == b2 || bead2 == b3)) ||
	(bead1 == b2 && (bead2 == b1 || bead2 == b3)) ||
	(bead1 == b3 && (bead2 == b1 || bead2 == b2))){
      notifyAngleRemoval(angle);
      bead1->removeAngle(angle);
      bead2->removeAngle(angle);
      angleList.erase(it);
      break;
    }
  }
}

void Bead::removeAllBonds(){
  for (auto const& bond : bondList){
    shared_ptr<Bead> bead1 = bond->getBead(0);
    shared_ptr<Bead> bead2 = bond->getBead(1);
    notifyBondRemoval(bond);
    if (bead1 != shared_from_this())
      bead1->removeBond(bond);
    else
      bead2->removeBond(bond);
  }
  bondList.clear();
}

void Bead::removeAllAngles(){
  for (auto const& angle : angleList){
    shared_ptr<Bead> bead1 = angle->getBead(0);
    shared_ptr<Bead> bead2 = angle->getBead(1);
    shared_ptr<Bead> bead3 = angle->getBead(2);
    notifyAngleRemoval(angle);
    if (bead1 == shared_from_this()){
      bead2->removeAngle(angle);
      bead3->removeAngle(angle);
    } else if (bead2 == shared_from_this()){
      bead1->removeAngle(angle);
      bead3->removeAngle(angle);
    } else {
      bead1->removeAngle(angle);
      bead2->removeAngle(angle);
    }
  }
  angleList.clear();  
}

// For handling bead, bond, and angle listeners
void Bead::addBeadListener(const shared_ptr<BeadListener>& l){
  auto it = std::find_if(beadListeners.begin(), beadListeners.end(),
			 [&](const weak_ptr<BeadListener>& p){
			   return p.lock() == l;});
  if (it == beadListeners.end()){
    beadListeners.push_back(l); 
  }
}

void Bead::removeBeadListener(const shared_ptr<BeadListener>& l){
  auto it = std::find_if(beadListeners.begin(), beadListeners.end(),
			 [&](const weak_ptr<BeadListener>& p){
			   return p.lock() == l;});
  if (it != beadListeners.end()){
    beadListeners.erase(it);
  }
}

void Bead::addBondListener(const shared_ptr<BondListener>& l){
  auto it = std::find_if(bondListeners.begin(), bondListeners.end(),
			 [&](const weak_ptr<BondListener>& p){
			   return p.lock() == l;});
  if (it == bondListeners.end()){
    bondListeners.push_back(l);
    // Register the listener to all existing bonds
    for (auto const& b: bondList){
      b->addListener(l);
    }
  }
}

void Bead::removeBondListener(const shared_ptr<BondListener>& l){
  auto it = std::find_if(bondListeners.begin(), bondListeners.end(),
			 [&](const weak_ptr<BondListener>& p){
			   return p.lock() == l;});
  if (it != bondListeners.end()){
    // De-register the listener from all existing bonds
    for (auto const& b : bondList){
      b->removeListener(it->lock());
    }
    bondListeners.erase(it);
  }
}

void Bead::addAngleListener(const shared_ptr<AngleListener>& l){
  auto it = std::find_if(angleListeners.begin(), angleListeners.end(),
			 [&](const weak_ptr<AngleListener>& p){
			   return p.lock() == l;});
  if (it == angleListeners.end()){
    angleListeners.push_back(l);
    // Register the listener to all existing angles
    for (auto const& a: angleList){
      a->addListener(l);
    }
  }
}

void Bead::removeAngleListener(const shared_ptr<AngleListener>& l){
  auto it = std::find_if(angleListeners.begin(), angleListeners.end(),
			 [&](const weak_ptr<AngleListener>& p){
			   return p.lock() == l;});
  if (it != angleListeners.end()){
    // De-register the listener from all existing angles
    for (auto const& a : angleList){
      a->removeListener(it->lock());
    }
    angleListeners.erase(it);
  }
}


void Bead::notifyBeadTypeChange(int oldType, int newType){
  for (auto const& l : beadListeners){
    l.lock()->beadTypeChanged(shared_from_this(), oldType, newType);
  }
}

void Bead::notifyBeadLabelChange(int oldLabel, int newLabel){
  for (auto const& l : beadListeners){
    l.lock()->beadLabelChanged(shared_from_this(), oldLabel, newLabel);
  }
}

void Bead::notifyBondCreation(const shared_ptr<Bond>& b){
  for (auto const& l : bondListeners){
    l.lock()->bondCreated(b);
    b->addListener(l.lock());
  }
}

void Bead::notifyBondRemoval(const shared_ptr<Bond>& b){
  for (auto const& l : bondListeners){
    b->removeListener(l.lock());
    l.lock()->bondRemoved(b);
  }
}

void Bead::notifyAngleCreation(const shared_ptr<Angle>& a){
  for (auto const& l : angleListeners){
    l.lock()->angleCreated(a);
    a->addListener(l.lock());
  }
}

void Bead::notifyAngleRemoval(const shared_ptr<Angle>& a){
  for (auto const& l : angleListeners){
    a->removeListener(l.lock());
    l.lock()->angleRemoved(a);
  }
}
