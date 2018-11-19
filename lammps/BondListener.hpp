/* BondListener.hpp
 *
 * An interface for handling bond events
 */

#ifndef BONDLISTENER_HPP
#define BONDLISTENER_HPP

#include <memory>
#include "SharedFromThis.hpp"

using std::shared_ptr;

class Bond; 

class BondListener :
  public inheritable_enable_shared_from_this<BondListener> {

public:
  virtual void bondCreated(const shared_ptr<Bond>& bond) = 0;
  virtual void bondRemoved(const shared_ptr<Bond>& bond) = 0;
  virtual void bondTypeChanged(const shared_ptr<Bond>& bond, 
			       int oldType, int newType) = 0;

};

#endif
