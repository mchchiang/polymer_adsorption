// LAMMPS.cpp

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <memory>
#include <map>
#include <string>
#include "Bead.hpp"
#include "Bond.hpp"
#include "Angle.hpp"
#include "Polymer.hpp"
#include "LAMMPS.hpp"
#include "BeadListener.hpp"
#include "BondListener.hpp"
#include "AngleListener.hpp"

using std::cout;
using std::endl;
using std::setw;
using std::setprecision;
using std::vector;
using std::string;
using std::map;
using std::shared_ptr;
using std::make_shared;
using std::ofstream;
using std::ifstream;
using std::istringstream;
using std::getline;

// Constructors
LAMMPS::LAMMPS() {}

LAMMPS::LAMMPS(double x, double y, double z) : lx {x}, ly {y}, lz {z} {}

// Accessor methods
shared_ptr<Polymer> LAMMPS::getPolymer(int id){
  return polymers[id];
}

shared_ptr<Bead> LAMMPS::getBead(int id){
  return beads[id];
}

void LAMMPS::setLx(double x){
  lx = x;
}

double LAMMPS::getLx(){
  return lx;
}

void LAMMPS::setLy(double y){
  ly = y;
}

double LAMMPS::getLy(){
  return ly;
}

void LAMMPS::setLz(double z){
  lz = z;
}

double LAMMPS::getLz(){
  return lz;
}

int LAMMPS::getNumOfBeads(){
  /*  int total {};
  for (auto const& p : polymers){
    total += p.second->getNumOfBeads();
  }
  total += beads.size();
  return total;*/
  return numOfBeads;
}

int LAMMPS::getNumOfBonds(){
  /*  int nbeads {}, total {};
  shared_ptr<Polymer> polymer;
  shared_ptr<Bead> bead;
  for (auto const& p : polymers){
	polymer = p.second;
	nbeads = polymer->getNumOfBeads();
	for (int i {}; i < nbeads; i++){
	  total += polymer->getBead(i)->getNumOfBonds();
	}
  }
  for (auto const& b : beads){
	bead = b.second;
	total += bead->getNumOfBonds();
  }
  return total/2;*/
  return numOfBonds;
}

int LAMMPS::getNumOfAngles(){
  /*  int nbeads {}, total {};
  shared_ptr<Polymer> polymer;
  shared_ptr<Bead> bead;
  for (auto const& p : polymers){
	polymer = p.second;
	nbeads = polymer->getNumOfBeads();
	for (int i {}; i < nbeads; i++){
	  total += polymer->getBead(i)->getNumOfAngles();
	}
  }
  for (auto const& b : beads){
	bead = b.second;
	total += bead->getNumOfAngles();
  }
  return total/3;*/
  return numOfAngles;
}

void LAMMPS::setTypesOfBeads(int type){
  typesOfBeads = type;
}

int LAMMPS::getTypesOfBeads(){
  //  return typesOfBeads;
  return beadTypeCount.size();
}

void LAMMPS::setTypesOfBonds(int type){
  typesOfBonds = type;
}

int LAMMPS::getTypesOfBonds(){
  //  return typesOfBonds;
  return bondTypeCount.size();
}

void LAMMPS::setTypesOfAngles(int type){
  typesOfAngles = type;
}

int LAMMPS::getTypesOfAngles(){
  //  return typesOfAngles;
  return angleTypeCount.size();
}

void LAMMPS::addBead(int id, shared_ptr<Bead> bead){
  beads[id] = bead;
}

void LAMMPS::removeBead(int id){
  if (beads.count(id) != 0){
    beads[id]->removeAllBonds();
    beads[id]->removeAllAngles();
    beads.erase(id);
  }
}

void LAMMPS::removeAllBeads(){
  for (auto const& b : beads){
    b.second->removeAllBonds();
    b.second->removeAllAngles();
  }
  beads.clear();
}

void LAMMPS::addPolymer(int id, shared_ptr<Polymer> polymer){
  polymers[id] = polymer;
}

void LAMMPS::removePolymer(int id){
  if (polymers.count(id) != 0){
    polymers[id]->removeAllBeads();
    polymers.erase(id);
  }
}

void LAMMPS::removeAllPolymers(){
  for (auto const& p : polymers){
    p.second->removeAllBeads();
  }
  polymers.clear();
}

void LAMMPS::clear(){
  removeAllPolymers();
  removeAllBeads();
}

bool LAMMPS::changeBeadID(int oldID, int newID){
  if (beads.count(oldID) != 0 && beads.count(newID) == 0){
    beads[newID] = beads[oldID];
    beads.erase(oldID);
    return true;
  }
  return false;
}

bool LAMMPS::changePolymerID(int oldID, int newID){
  if (polymers.count(oldID) != 0 && polymers.count(newID) == 0){
    polymers[newID] = polymers[oldID];
    polymers.erase(oldID);
    return true;
  }
  return false;
}

bool LAMMPS::importData(string inFile, string mapFile){
  // Clear currently stored data
  clear();

  cout << "Start reading LAMMPS file ..." << endl;
  
  ifstream reader;
  reader.open(inFile);
  
  int numOfBeads {};
  int numOfBonds {};
  int numOfAngles {};
  string line {};

  if (!readHeader(reader, numOfBeads, numOfBonds, numOfAngles)){
    reader.close();
    cout << "Data file must specify number and types of "
	 << "atoms, bonds, and angles." << endl;
    return false;
  }
  
  if (!readBoxSize(reader)){
    reader.close();
    cout << "Data file must specify box size in all three dimensions." << endl;
    return false;
  }
  
  bool positionOK {false};
  bool velocityOK {false};
  bool bondOK {false};
  bool angleOK {false};

  map< int, shared_ptr<Bead> > beadIndexMap {};
  
  if (numOfBonds == 0) bondOK = true;
  if (numOfAngles == 0) angleOK = true;

  // Create all the beads
  for (int i {1}; i <= numOfBeads; i++){
    beadIndexMap[i] = make_shared<Bead>();
  }

  // Read position, velocity, bond, and angle data
  while (!reader.eof() && (!positionOK || !velocityOK || !bondOK || !angleOK)){
    getline(reader, line);
	
    // Read atoms' positions
    if (line.find("Atoms") != string::npos && !positionOK){
      positionOK = readPosition(reader, numOfBeads, beadIndexMap);

      // Read atoms' velocities
    } else if (line.find("Velocities") != string::npos && !velocityOK){
      velocityOK = readVelocity(reader, numOfBeads, beadIndexMap);

      // Read bonds
    } else if (line.find("Bonds") != string::npos && !bondOK){
      bondOK = readBond(reader, numOfBonds, beadIndexMap);

      // Read angles
    } else if (line.find("Angles") != string::npos && !angleOK){
      angleOK = readAngle(reader, numOfAngles, beadIndexMap);
    }
  }

  if (!positionOK || !velocityOK || !bondOK || !angleOK){
    reader.close();
    cout << "Problem with reading position, "
	 << "velocity, bond, or angle data" << endl;
    return false;
  }
  
  // Read the mapping file
  if (mapFile.compare("") != 0){

    cout << "Reading mapping file" << endl;

    ifstream mapReader;
    mapReader.open(mapFile);
    readInputMap(mapReader, beadIndexMap);
    mapReader.close();

  }

  // Add remaining beads to container
  for (auto const& b : beadIndexMap){
    beads[b.first] = b.second;
  }
  
  reader.close();
  
  cout << "Finish reading LAMMPS file" << endl;

  return true;
}

bool LAMMPS::readHeader(ifstream& reader, int& numOfBeads,
			int& numOfBonds, int& numOfAngles){

  cout << "Reading simulation info ..." << endl;

  bool readNumOfBeads {false};
  bool readTypesOfBeads {false};
  bool readNumOfBonds {false};
  bool readTypesOfBonds {false};
  bool readNumOfAngles {false};
  bool readTypesOfAngles {false};

  string line {};
  istringstream ss;
  
  while (!readNumOfBeads || !readTypesOfBeads ||
	 !readNumOfBonds || !readTypesOfBonds ||
	 !readNumOfAngles || !readTypesOfAngles){
    if (reader.eof()) return false;
    getline(reader, line);
    ss.clear();
    ss.str(line);
    
    if (line.find("atoms") != string::npos){
      ss >> numOfBeads;
      readNumOfBeads = true;
    } else if (line.find("atom types") != string::npos){
      ss >> typesOfBeads;
      readTypesOfBeads = true;
    } else if (line.find("bonds") != string::npos){
      ss >> numOfBonds;
      readNumOfBonds = true;
    } else if (line.find("bond types") != string::npos){
      ss >> typesOfBonds;
      readTypesOfBonds = true;
    } else if (line.find("angles") != string::npos){
      ss >> numOfAngles;
      readNumOfAngles = true;
    } else if (line.find("angle types") != string::npos){
      ss >> typesOfAngles;
      readTypesOfAngles = true;
    }
  }

  cout << numOfBeads << " atoms" << endl;
  cout << typesOfBeads << " atom types" << endl;
  cout << numOfBonds << " bonds" << endl;
  cout << typesOfBonds << " bond types" << endl;
  cout << numOfAngles << " angles" << endl;
  cout << typesOfAngles << " angle types" << endl;
  
  return true;
}

bool LAMMPS::readBoxSize(ifstream& reader){
  
  cout << "Reading box size ... " << endl;
  
  bool readLx {false};
  bool readLy {false};
  bool readLz {false};
  string line {};
  istringstream ss;
  
  double lleft {}, lright {};

  if (reader.eof()) return false;
  
  while (!readLx || !readLy || !readLz){
    getline(reader, line);
    ss.clear();
    ss.str(line);

    if (line.find("xlo xhi") != string::npos){
      ss >> lleft >> lright;
      lx = lright - lleft;
      readLx = true;
    } else if (line.find("ylo yhi") != string::npos){
      ss >> lleft >> lright;
      ly = lright - lleft;
      readLy = true;
    } else if (line.find("zlo zhi") != string::npos){
      ss >> lleft >> lright;
      lz = lright - lleft;
      readLz = true;
    }

    if (reader.eof()) return false;
  }
  
  cout << "Lx = " << lx << endl;
  cout << "Ly = " << ly << endl;
  cout << "Lz = " << lz << endl;
  
  return true;
}

bool LAMMPS::readPosition(ifstream& reader, int& numOfBeads, 
			  map< int, shared_ptr<Bead> >& beadIndexMap){
  int count {}, index {}, label {}, type {};
  double x {}, y {}, z {};
  int nx {}, ny {}, nz {};
  shared_ptr<Bead> bead {};
  string line {};
  istringstream ss;

  cout << "Reading atoms' positions ..." << endl;
  
  while (count < numOfBeads && !reader.eof()){
    getline(reader, line);
    if (line.length() > 0 && line[0] != '#'){
      ss.clear();
      ss.str(line);
      ss >> index >> label >> type >> x >> y >> z >> nx >> ny >> nz;
      bead = beadIndexMap[index];
      bead->setLabel(label);
      bead->setType(type);
      bead->setPosition(0, x);
      bead->setPosition(1, y);
      bead->setPosition(2, z);
      bead->setBoundaryCount(0, nx);
      bead->setBoundaryCount(1, ny);
      bead->setBoundaryCount(2, nz);		  
      count++;
    }
  }

  if (count == numOfBeads){
    return true;
  }
  return false;
}

bool LAMMPS::readVelocity(ifstream& reader, int& numOfBeads,
			  map< int, shared_ptr<Bead> >& beadIndexMap){
  int count {}, index {};
  double vx {}, vy {}, vz {};
  shared_ptr<Bead> bead {};
  string line {};
  istringstream ss;
  
  cout << "Reading atoms' velocities ..." << endl;
  
  while (count < numOfBeads && !reader.eof()){
    getline(reader, line);
    if (line.length() > 0 && line[0] != '#'){
      ss.clear();
      ss.str(line);
      ss >> index >> vx >> vy >> vz;
      bead = beadIndexMap[index];
      bead->setVelocity(0, vx);
      bead->setVelocity(1, vy);
      bead->setVelocity(2, vz);
      count++;
    }
  }
  
  if (count == numOfBeads){
    return true;
  }
  return false;
}

bool LAMMPS::readBond(ifstream& reader, int& numOfBonds,
		      map< int, shared_ptr<Bead> >& beadIndexMap){
  int count {}, index {}, type {};
  int bead1Index {}, bead2Index {};
  shared_ptr<Bead> bead1 {}, bead2 {};
  string line {};
  istringstream ss;

  cout << "Reading bonds ..." << endl;
  
  while (count < numOfBonds && !reader.eof()){
    getline(reader, line);
    if (line.length() > 0 && line[0] != '#'){
      ss.clear();
      ss.str(line);
      ss >> index >> type >> bead1Index >> bead2Index;
      bead1 = beadIndexMap[bead1Index];
      bead2 = beadIndexMap[bead2Index];
      bead1->addBondWith(type, bead2);
      count++;
    }
  }
  
  if (count == numOfBonds){
    return true;
  }
  return false;
}

bool LAMMPS::readAngle(ifstream& reader, int& numOfAngles,
		       map< int, shared_ptr<Bead> >& beadIndexMap){
  int count {}, index {}, type {};
  int bead1Index {}, bead2Index {}, bead3Index {};
  shared_ptr<Bead> bead1 {}, bead2 {}, bead3 {};
  string line {};
  istringstream ss;

  cout << "Reading angles ..." << endl;
  
  while (count < numOfAngles && !reader.eof()){
    getline(reader, line);
    if (line.length() > 0 && line[0] != '#'){
      ss.clear();
      ss.str(line);
      ss >> index >> type >> bead1Index >> bead2Index >> bead3Index;
      bead1 = beadIndexMap[bead1Index];
      bead2 = beadIndexMap[bead2Index];
      bead3 = beadIndexMap[bead3Index];
      bead1->addAngleWith(type, bead2, bead3);
      count++;
    }
  }
  
  if (count == numOfAngles){
    return true;
  }
  return false;
}

void LAMMPS::readInputMap(ifstream& reader, 
			  map< int, shared_ptr<Bead> >& beadIndexMap){
  istringstream ss;
  string line {};
  
  bool readingPolymer {false};
  bool readingBead {false};
  
  int size {}, key {}, startBead {}, endBead {};
  shared_ptr<Polymer> polymer {};
  
  while (!reader.eof()){
    getline(reader, line);
    ss.clear();
    ss.str(line);
	
    if (line.length() > 0 && line[0] == '#'){
      if (line.
	  compare("# Polymer mapping - polymer key, start bead, end bead")
	  == 0){
	readingPolymer = true;
	readingBead = false;
      } else if (line.
		 compare("# Bead mapping - bead key, bead index") == 0){
	readingBead = true;
	readingPolymer = false;
      }
    } else if (line.length() > 0) {
      if (readingPolymer){
	ss >> key >> startBead >> endBead;
	size = endBead - startBead + 1;
	polymer = make_shared<Polymer>(size, 1, 1, 1, false);
	for (int i {startBead}; i <= endBead; i++){
	  polymer->addBead(beadIndexMap[i]);
	  beadIndexMap.erase(i);
	}
	polymers[key] = polymer;
      } else if (readingBead){
	ss >> key >> startBead;
	beads[key] = beadIndexMap[startBead];
	beadIndexMap.erase(startBead);
      }
    }
  }
}

bool LAMMPS::exportData(string outFile, string mapFile){
  int beadIndexCount {1};
  int bondIndexCount {1};
  int angleIndexCount {1};

  map< shared_ptr<Bead>, int> beadIndexMap {};
  map< shared_ptr<Bond>, int> bondIndexMap {};
  map< shared_ptr<Angle>, int> angleIndexMap {};  

  cout << "Start writing LAMMPS file ..." << endl;

  // Writer for the LAMMPS input file
  ofstream writer;
  writer.open(outFile);

  stringstream headerWriter;
  stringstream positionWriter;
  stringstream velocityWriter;
  stringstream bondWriter;
  stringstream angleWriter;

  positionWriter << "\nAtoms\n" << endl;
  positionWriter.unsetf(std::ios_base::floatfield);
  // positionWriter << std::defaultfloat;
  velocityWriter << "\nVelocities\n" << endl;
  velocityWriter.unsetf(std::ios_base::floatfield);
  // velocityWriter << std::defaultfloat;
  bondWriter << "\nBonds\n" << endl;
  bondWriter.unsetf(std::ios_base::floatfield);
  // bondWriter << std::defaultfloat;
  angleWriter << "\nAngles\n" << endl;
  angleWriter.unsetf(std::ios_base::floatfield);
  // angleWriter << std::defaultfloat;

  // Writer for the mapping between indices and polymers/beads
  ofstream mapWriter;
  mapWriter.open(mapFile);
  mapWriter << "# Mapping file between bead indices and polymers" << endl;
  mapWriter << "# Polymer mapping - polymer key, start bead, end bead" << endl;

  // Write positions and velocities
  for (auto const& p : polymers){
    mapWriter << p.first << " " << beadIndexCount << " ";
    for (auto const& b : p.second->getBeads()){
      writePositionAndVelocity(b, beadIndexMap, 
			       positionWriter, velocityWriter, 
			       beadIndexCount);
    }
    mapWriter << beadIndexCount-1 << endl;
  }

  mapWriter << endl;
  mapWriter << "# Bead mapping - bead key, bead index" << endl;

  for (auto const& b : beads){
    mapWriter << b.first << " " << beadIndexCount << endl;
    writePositionAndVelocity(b.second, beadIndexMap, 
			     positionWriter, velocityWriter, 
			     beadIndexCount);
  }
  
  // Write bonds and angles
  for (auto const& p : polymers){
    for (auto const& b : p.second->getBeads()){
      writeBondAndAngle(b, beadIndexMap, bondIndexMap, angleIndexMap,
			bondWriter, angleWriter,
			bondIndexCount, angleIndexCount);
    }
  }
 
  for (auto const& b : beads){
    writeBondAndAngle(b.second, beadIndexMap, bondIndexMap, angleIndexMap,
		      bondWriter, angleWriter,
		      bondIndexCount, angleIndexCount);
  }

  beadIndexCount--;
  bondIndexCount--;
  angleIndexCount--;

  writeHeader(headerWriter,
	      beadIndexCount, bondIndexCount, angleIndexCount);
  
  headerWriter << endl;
  positionWriter << endl;
  velocityWriter << endl;
  bondWriter << endl;
  angleWriter << endl;

  writer << headerWriter.str();
  writer << positionWriter.str();
  writer << velocityWriter.str();
  if (bondIndexCount > 0)
    writer << bondWriter.str();
  if (angleIndexCount > 0)
    writer << angleWriter.str();
  
  writer.close();
  
  cout << "Finish writing LAMMPS file" << endl;
  
  return true;
}

void LAMMPS::writePositionAndVelocity(const shared_ptr<Bead>& bead,
				      map<shared_ptr<Bead>, int>& beadIndexMap,
				      stringstream& positionWriter,
				      stringstream& velocityWriter,
				      int& beadIndexCount){
  if (beadIndexMap[bead] == 0){ // Bead is not in the index map
    writePosition(positionWriter, bead, beadIndexCount);
    writeVelocity(velocityWriter, bead, beadIndexCount);
    beadIndexMap[bead] = beadIndexCount;
    beadIndexCount++;
  }
}

void LAMMPS::writeBondAndAngle(const shared_ptr<Bead>& bead,
			       map<shared_ptr<Bead>,int>& beadIndexMap,
			       map<shared_ptr<Bond>,int>& bondIndexMap,
			       map<shared_ptr<Angle>,int>& angleIndexMap,
			       stringstream& bondWriter,
			       stringstream& angleWriter,
			       int& bondIndexCount, int& angleIndexCount){
  vector<shared_ptr<Bond> > bondList = bead->getBonds();
  vector<shared_ptr<Angle> > angleList = bead->getAngles();
  
  for (auto const& bond : bondList){
    if (bondIndexMap[bond] == 0){ // Bond is not in the index map
      bondIndexMap[bond] = bondIndexCount;
      int type = bond->getType();
      int bead1Index = beadIndexMap[bond->getBead(0)];
      int bead2Index = beadIndexMap[bond->getBead(1)];
      writeBond(bondWriter, bondIndexCount, type, bead1Index, bead2Index);
      bondIndexCount++;
    }
  }
  
  for (auto const& angle : angleList){
    if (angleIndexMap[angle] == 0){ // Angle is not in the index map
      angleIndexMap[angle] = angleIndexCount;
      int type = angle->getType();
      int bead1Index = beadIndexMap[angle->getBead(0)];
      int bead2Index = beadIndexMap[angle->getBead(1)];
      int bead3Index = beadIndexMap[angle->getBead(2)];
      writeAngle(angleWriter, angleIndexCount, type, 
		 bead1Index, bead2Index, bead3Index);
      angleIndexCount++;
    }
  }
} 

void LAMMPS::writeHeader(stringstream& writer,
			 int numOfBeads, int numOfBonds, int numOfAngles){
  string header {
    "LAMMPS data file from restart file: timestep = 0,\tprocs = 1"};
  
  writer << header << endl;
  writer << endl;
  writer << numOfBeads << " atoms " << endl;
  writer << numOfBonds << " bonds " << endl;
  writer << numOfAngles << " angles " << endl;
  writer << "\n";
  writer << typesOfBeads << " atom types " << endl;
  writer << typesOfBonds << " bond types " << endl;
  writer << typesOfAngles << " angle types " << endl;
  writer << "\n";
  writer << -lx/2.0 << " " << (lx-lx/2.0) << " xlo xhi" << endl;
  writer << -ly/2.0 << " " << (ly-ly/2.0) << " ylo yhi" << endl;
  writer << -lz/2.0 << " " << (lz-lz/2.0) << " zlo zhi" << endl;
  
  writer << "\nMasses\n" << endl;
  for (int i {1}; i <= typesOfBeads; i++){
    writer << i << " " << 1 << endl;
  }
}

void LAMMPS::writePosition(stringstream& writer,
			   const shared_ptr<Bead>& bead, int beadIndex){
  writer << beadIndex << " "
         << bead->getLabel() << " "
         << bead->getType() << " ";
  writer << std::scientific;
  writer << setprecision(preci) << bead->getPosition(0) << " "
         << setprecision(preci) << bead->getPosition(1) << " "
         << setprecision(preci) << bead->getPosition(2) << " ";
  writer.unsetf(std::ios_base::floatfield);
  // writer << std::defaultfloat;
  writer << bead->getBoundaryCount(0) << " "
         << bead->getBoundaryCount(1) << " "
         << bead->getBoundaryCount(2) << endl;
}

void LAMMPS::writeVelocity(stringstream& writer,
			   const shared_ptr<Bead>& bead, int beadIndex){
  writer << beadIndex << " ";
  writer << std::scientific;
  writer << setprecision(preci) << bead->getVelocity(0) << " "
         << setprecision(preci) << bead->getVelocity(1) << " "
         << setprecision(preci) << bead->getVelocity(2) << endl;
  writer.unsetf(std::ios_base::floatfield);
  // writer << std::defaultfloat;
}

void LAMMPS::writeBond(stringstream& writer, int bondIndex, int bondType,
		       int bead1Index, int bead2Index){
  writer << bondIndex << " "
         << bondType << " "
         << bead1Index << " "
         << bead2Index << endl;
}

void LAMMPS::writeAngle(stringstream& writer, int angleIndex, int angleType,
			int bead1Index, int bead2Index, int bead3Index){
  writer << angleIndex << " "
         << angleType  << " "
         << bead1Index << " "
         << bead2Index << " "
         << bead3Index << endl;
}

// For handling bead, bond, and angle events
void LAMMPS::beadTypeChanged(const shared_ptr<Bead>& bead, 
			     int oldType, int newType){
  beadTypeCount[oldType]--;
  beadTypeCount[newType]++;
  if (beadTypeCount[oldType] <= 0){
    beadTypeCount.erase(oldType);
  }
}

void LAMMPS::beadLabelChanged(const shared_ptr<Bead>& bead,
			      int oldLabel, int newLabel){
  // LAMMPS don't need to no bead label changes
}

void LAMMPS::bondCreated(const shared_ptr<Bond>& bond){
  numOfBonds++;
  bondTypeCount[bond->getType()]++;
}
  
void LAMMPS::bondRemoved(const shared_ptr<Bond>& bond){
  numOfBonds--;
  int type {bond->getType()};
  bondTypeCount[type]--;
  if (bondTypeCount[type] <= 0){
    bondTypeCount.erase(type);
  }
}

void LAMMPS::bondTypeChanged(const shared_ptr<Bond>& bond, 
			     int oldType, int newType){
  bondTypeCount[oldType]--;
  bondTypeCount[newType]++;
  if (bondTypeCount[oldType] <= 0){
    bondTypeCount.erase(oldType);
  }
}

void LAMMPS::angleCreated(const shared_ptr<Angle>& angle){
  numOfAngles++;
  angleTypeCount[angle->getType()]++;
}

void LAMMPS::angleRemoved(const shared_ptr<Angle>& angle){
  numOfAngles--;
  int type {angle->getType()};
  angleTypeCount[type]--;
  if (angleTypeCount[type] <= 0){
    angleTypeCount.erase(type);
  }
}

void LAMMPS::angleTypeChanged(const shared_ptr<Angle>& angle,
			      int oldType, int newType){
  angleTypeCount[oldType]--;
  angleTypeCount[newType]++;
  if (angleTypeCount[oldType] <= 0){
    angleTypeCount.erase(oldType);
  }
}

void LAMMPS::beadCreated(const shared_ptr<Bead>& bead, int beadType){
  bead->addBeadListener(BeadListener::shared_from_this());
  bead->addBondListener(BondListener::shared_from_this());
  bead->addAngleListener(AngleListener::shared_from_this());
  numOfBeads++;
  beadTypeCount[beadType]++;
}

void LAMMPS::polymerCreated(const shared_ptr<Polymer>& polymer, int nBeads,
			    int beadType, int bondType, int angleType){
  polymer->addBeadListener(BeadListener::shared_from_this());
  polymer->addBondListener(BondListener::shared_from_this());
  polymer->addAngleListener(AngleListener::shared_from_this());
  numOfBeads += nBeads;
  numOfBonds += nBeads-1;
  numOfAngles += nBeads-2;
  beadTypeCount[beadType] += nBeads;
  bondTypeCount[bondType] += nBeads-1;
  angleTypeCount[angleType] += nBeads-2;
}

// For adding polymers and beads to the system
shared_ptr<Polymer> LAMMPS::createPolymer(int id, int nBeads,
					  int beadType, int bondType,
					  int angleType){
  removePolymer(id);
  shared_ptr<Polymer> polymer = make_shared<Polymer>(nBeads, beadType,
						     bondType, angleType);
  polymerCreated(polymer, nBeads, beadType, bondType, angleType);
  addPolymer(id, polymer);
  return polymer;
}

shared_ptr<Polymer> 
LAMMPS::createRandomWalkPolymer(int id, int nBeads, int beadType,
				int bondType, int angleType,
				double x0, double y0, double z0, 
				double rx, double ry, double rz){
  removePolymer(id);
  shared_ptr<Polymer> polymer = 
    Polymer::createRandomWalkPolymer(nBeads, beadType, bondType, angleType,
				     x0, y0, z0, rx, ry, rz);
  polymerCreated(polymer, nBeads, beadType, bondType, angleType);
  addPolymer(id, polymer);
  return polymer;
}

shared_ptr<Bead> LAMMPS::createBead(int id, int beadType, int beadLabel){
  removeBead(id);
  shared_ptr<Bead> bead = make_shared<Bead>(0,0,0,0,0,0,0,0,0,
					    beadType, beadLabel);
  beadCreated(bead, beadType);
  addBead(id, bead);
  return bead;
}

shared_ptr<Bead> LAMMPS::createBead(int id, double x, double y, double z,
				    double vx, double vy, double vz,
				    double nx, double ny, double nz,
				    int beadType, int beadLabel){
  removeBead(id);
  shared_ptr<Bead> bead = make_shared<Bead>(x,y,z,vx,vy,vz,nx,ny,nz,
					    beadType, beadLabel);
  beadCreated(bead, beadType);
  addBead(id, bead);
  return bead;
}
