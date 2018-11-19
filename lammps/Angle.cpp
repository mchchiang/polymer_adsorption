// Angle.cpp

#include <memory>
#include <algorithm>
#include "Bead.hpp"
#include "Angle.hpp"
#include "AngleListener.hpp"

using std::weak_ptr;
using std::shared_ptr;

// Constructor
Angle::Angle(int t, const shared_ptr<Bead>& b1, 
	     const shared_ptr<Bead>& b2, const shared_ptr<Bead>& b3) :
  type {t}, beads {b1,b2,b3} {}

// Accessor methods
shared_ptr<Bead> Angle::getBead(int id){
  return beads[id].lock();
}

void Angle::setType(int t){
  int oldType {type};
  type = t;
  if (oldType != type){
    notifyTypeChange(oldType, type);
  }
}

int Angle::getType(){
  return type;
}

// For handling listeners
void Angle::addListener(const shared_ptr<AngleListener>& l){
  auto it = std::find_if(listeners.begin(), listeners.end(),
			 [&](const weak_ptr<AngleListener>& p){
			   return p.lock() == l;});
  if (it == listeners.end()){
    listeners.push_back(l);
  }
}

void Angle::removeListener(const shared_ptr<AngleListener>& l){
  auto it = std::find_if(listeners.begin(), listeners.end(),
			 [&](const weak_ptr<AngleListener>& p){
			   return p.lock() == l;});
  if (it != listeners.end()){
    listeners.erase(it);
  }
}

void Angle::notifyTypeChange(int oldType, int newType){
  for (auto const& l : listeners){
    l.lock()->angleTypeChanged(shared_from_this(), oldType, newType);
  }
}
