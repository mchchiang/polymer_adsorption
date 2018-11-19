// Bond.cpp

#include <memory>
#include <algorithm>
#include "Bead.hpp"
#include "Bond.hpp"
#include "BondListener.hpp"

using std::weak_ptr;
using std::shared_ptr;

// Constructor
Bond::Bond(int t, const shared_ptr<Bead>& b1, const shared_ptr<Bead>& b2) :
  type {t}, beads {b1,b2} {}

// Accessor methods
shared_ptr<Bead> Bond::getBead(int id){
  return beads[id].lock();
}

void Bond::setType(int t){
  int oldType {type};
  type = t;
  if (oldType != type){
    notifyTypeChange(oldType, type);
  }
}

int Bond::getType(){
  return type;
}

// For handling listeners
void Bond::addListener(const shared_ptr<BondListener>& l){
  auto it = std::find_if(listeners.begin(), listeners.end(),
			 [&](const weak_ptr<BondListener>& p){
			   return p.lock() == l;});
  if (it == listeners.end()){
    listeners.push_back(l);
  }
}

void Bond::removeListener(const shared_ptr<BondListener>& l){
  auto it = std::find_if(listeners.begin(), listeners.end(),
			 [&](const weak_ptr<BondListener>& p){
			   return p.lock() == l;});
  if (it != listeners.end()){
    listeners.erase(it);
  }
}

void Bond::notifyTypeChange(int oldType, int newType){
  for (auto const& l : listeners){
    l.lock()->bondTypeChanged(shared_from_this(), oldType, newType);
  }
}
