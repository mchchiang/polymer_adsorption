/* AngleListener.hpp
 *
 * An interface for handling bond events
 */

#ifndef ANGLELISTENER_HPP
#define ANGLELISTENER_HPP

#include <memory>
#include "SharedFromThis.hpp"

using std::shared_ptr;

class Angle; 

class AngleListener : 
  public inheritable_enable_shared_from_this<AngleListener> {

public:
  virtual void angleCreated(const shared_ptr<Angle>& angle) = 0;
  virtual void angleRemoved(const shared_ptr<Angle>& angle) = 0;
  virtual void angleTypeChanged(const shared_ptr<Angle>& angle, 
			       int oldType, int newType) = 0;

};

#endif
