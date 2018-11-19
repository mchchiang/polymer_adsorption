// PositionReader.cpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "PositionReader.hpp"

using std::cout;
using std::endl;
using std::ifstream;
using std::istringstream;
using std::vector;
using std::string;


// Destructor
PositionReader::~PositionReader(){
  deleteData();
}

void PositionReader::initData(){
  vector<double> zeroDouble (3, 0.0);
  vector<int> zeroInt (3, 0);
  position = new vector<vector<double> >(numOfBeads, zeroDouble);
  boundaryCount = new vector<vector<int> >(numOfBeads, zeroInt);
  type = new vector<int>(numOfBeads, 0);
}

void PositionReader::deleteData(){
  if (!position) delete position;
  if (!boundaryCount) delete boundaryCount;
  if (!type) delete type;
}

bool PositionReader::open(string posFile, int nBeads,
			  double lx, double ly, double lz, int tInc){
  if (fileOpen){
    cout << "Close current file before opening another file" << endl;
    return false;
  }
  reader.open(posFile);
  if (!reader){
    cout << "Problem with opening position file!" << endl;
    return false;
  }
  fileOpen = true;
  numOfBeads = nBeads;
  boxSize = {lx, ly, lz};
  timeInc = tInc;
  time = -timeInc;
  initData();
  return true;
}

void PositionReader::close(){
  if (reader.is_open()){
    reader.close();
  }
  deleteData();
  fileOpen = false;
}

double PositionReader::getPosition(int beadIndex, int comp) const{
  return (*position)[beadIndex][comp];
}

double PositionReader::getUnwrappedPosition(int beadIndex, int comp) const{
  return (*position)[beadIndex][comp] + 
    boxSize[comp] * (*boundaryCount)[beadIndex][comp];
}

int PositionReader::getType(int beadIndex) const{
  return (*type)[beadIndex];
}

bool PositionReader::nextFrame(){
  if (!fileOpen || reader.eof()){
    return false;
  }
  string line, sym;
  double x, y, z;
  int ix, iy, iz, t;
  istringstream iss;
  const int headerLines {2};
  
  // Ignore header lines
  for (int i {}; i < headerLines; i++){
    if (!getline(reader, line)) return false;
  }

  // Read bead position data
  for (int i {}; i < numOfBeads; i++){
    if (!getline(reader, line)) return false;
    iss.clear();
    iss.str(line);
    iss >> sym >> x >> y >> z >> ix >> iy >> iz >> t;
    (*position)[i][0] = x;
    (*position)[i][1] = y;
    (*position)[i][2] = z;
    (*boundaryCount)[i][0] = ix;
    (*boundaryCount)[i][1] = iy;
    (*boundaryCount)[i][2] = iz;
    (*type)[i] = t;
  }

  time += timeInc;
  return true;
}

long PositionReader::getTime() const{
  return time;
}

bool PositionReader::isOpen(){
  return fileOpen;
}
