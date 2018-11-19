/* BeadListener.hpp
 *
 * An interface for handling bond events
 */

#ifndef BEADLISTENER_HPP
#define BEADLISTENER_HPP

#include <memory>
#include "SharedFromThis.hpp"

using std::shared_ptr;

class Bead;
class Bond; 
class Angle;

class BeadListener : 
  public inheritable_enable_shared_from_this<BeadListener> {

public:
  virtual void beadTypeChanged(const shared_ptr<Bead>& bead, 
			       int oldType, int newType) = 0;
  virtual void beadLabelChanged(const shared_ptr<Bead>& bead,
				int oldLabel, int newLabel) = 0;
};

#endif
