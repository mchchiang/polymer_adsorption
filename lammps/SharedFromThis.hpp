/* 
 * SharedFromThis.hpp
 * To allow multiple inheritance of objects inheriting shared_from_this.
 */

#ifndef SHAREDFROMTHIS_HPP
#define SHAREDFROMTHIS_HPP

#include <memory>

// A common base class
class MultiInheritableEnableSharedFromThis : 
  public std::enable_shared_from_this<MultiInheritableEnableSharedFromThis> {
public:
  virtual ~MultiInheritableEnableSharedFromThis() {}
};

template <class T>
class inheritable_enable_shared_from_this : 
  virtual public MultiInheritableEnableSharedFromThis {
public:
  std::shared_ptr<T> shared_from_this() {
    return std::dynamic_pointer_cast<T>(MultiInheritableEnableSharedFromThis::
					shared_from_this());
  }

  template <class Down>
  std::shared_ptr<Down> downcasted_shared_from_this(){
    return 
      std::dynamic_pointer_cast<Down>(MultiInheritableEnableSharedFromThis::
				      shared_from_this());
  }
};

#endif
