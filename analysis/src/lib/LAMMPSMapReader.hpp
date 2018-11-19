/* MapReader.hpp
 * A code which reads the lammpsmap file
 */

#ifndef LAMMPSMAPREADER_HPP
#define LAMMPSMAPREADER_HPP

#include <fstream>
#include <map>
#include <vector>
#include <string>

class LAMMPSMapReader {
  
private:
  bool fileOpen {false};
  std::ifstream reader;
  std::map<int,std::vector<int> > polymerMap;
  std::map<int,int> beadMap;

  void initData();

public:
  // Load a new map file
  bool open(std::string mapFile);
  void close();
  bool isOpen();
  
  std::map<int,int> getBeadMap();
  std::map<int,std::vector<int> > getPolymerMap();
  int getPolymer(int index);
  int getBead(int index);
  bool isPolymer(int index);
  bool isBead(int index);
  
};

#endif
