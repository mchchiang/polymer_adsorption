// LAMMPSMapReader.cpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>

#include "LAMMPSMapReader.hpp"

using std::cout;
using std::endl;
using std::istringstream;
using std::vector;
using std::map;
using std::string;
using LMReader = LAMMPSMapReader;

void LMReader::initData(){
  // Read the map file content
  istringstream iss;
  string line;
  int key, startBead, endBead;
  
  bool readingPolymer {false};
  bool readingBead {false};
  
  while (!reader.eof()){
    getline(reader, line);
    if (line.length() > 0 && line[0] == '#'){
      if (line.compare("# Polymer mapping - polymer key, start bead, end bead")
	  == 0){
	readingPolymer = true;
	readingBead = false;
      } else if (line.compare("# Bead mapping - bead key, bead index") == 0){
	readingBead = true;
	readingPolymer = false;
      }
    } else if (line.length() > 0){
      iss.clear();
      iss.str(line);
      if (readingPolymer){
	iss >> key >> startBead >> endBead;
	polymerMap[key] = {startBead-1, endBead-1};
      } else if (readingBead){
	iss >> key >> startBead;
	beadMap[key] = startBead-1;
      }
    }
  }
}

bool LMReader::open(string mapFile){
  if (fileOpen){
    cout << "Close current file before opening another file" << endl;
    return false;
  }
  reader.open(mapFile);
  if (!reader){
    cout << "Problem with opening map file!" << endl;
    return false;
  }
  fileOpen = true;
  initData();
  return true;
}

void LMReader::close(){
  reader.close();
}

bool LMReader::isOpen(){
  return fileOpen;
}

map<int,vector<int> > LMReader::getPolymerMap(){
  map<int,vector<int> > copy {polymerMap};
  return copy;
}

map<int,int> LMReader::getBeadMap(){
  map<int,int> copy {beadMap};
  return copy;
}

int LMReader::getPolymer(int index){
  for (auto const& p : polymerMap){
    if (index >= p.second[0] && index <= p.second[1]){
      return p.first;
    }
  }
  return -1;
}

int LMReader::getBead(int index){
  for (auto const& b : beadMap){
    if (index == b.second){
      return b.first;
    }
  }
  return -1;
}

bool LMReader::isPolymer(int index){
  for (auto const& p : polymerMap){
    if (index >= p.second[0] && index <= p.second[1]){
      return true;
    }
  }
  return false;
}

bool LMReader::isBead(int index){
  for (auto const& b : beadMap){
    if (index == b.second){
      return true;
    }
  }
  return false;
}
