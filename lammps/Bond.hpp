/* Bond.hpp
 *
 * A class for storing the info of two beads that are bonded
 */

#ifndef BOND_HPP
#define BOND_HPP

#include <vector>
#include <memory>
#include "BondListener.hpp"

using std::vector;
using std::weak_ptr;
using std::shared_ptr;

class Bead; // Forward declaration to break cyclic references

class Bond : public std::enable_shared_from_this<Bond> {
  
private:
  const int numOfBeads {2};
  int type;
  weak_ptr<Bead> beads[2];
  vector<weak_ptr<BondListener> > listeners {};

  void notifyTypeChange(int oldType, int newType);
  
public:
  // Constructor
  Bond(int type, const shared_ptr<Bead>& bead1, const shared_ptr<Bead>& bead2);
  
  // Accessor methods
  shared_ptr<Bead> getBead(int id);
  void setType(int type);
  int getType();

  // For handling listeners
  void addListener(const shared_ptr<BondListener>& listener);
  void removeListener(const shared_ptr<BondListener>& listener);
  
};

#endif
