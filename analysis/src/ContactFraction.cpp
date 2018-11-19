/* ContactFraction.cpp
 * A program that reads the lammpstrj file and determine
 * the fraction of beads in contact with the wall
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>  // for std::fill
#include "PositionReader.hpp"

using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::vector;

int main(int argc, char* argv[]){
  if (argc != 12){
    cout << "Usage: ContactFraction [numOfBeads] [typesOfBeads] "
	 << " [lx] [ly] [lz] [wallDist] "
	 << " [startTime] [endTime] [timeInc] [posFile] [outFile]" << endl;
    return 1;
  }

  int argi {};
  int numOfBeads {stoi(string(argv[++argi]), nullptr, 10)};
  int typesOfBeads {stoi(string(argv[++argi]), nullptr, 10)};
  double lx {stod(string(argv[++argi]), nullptr)};
  double ly {stod(string(argv[++argi]), nullptr)};
  double lz {stod(string(argv[++argi]), nullptr)};
  double wallDist {stod(string(argv[++argi]), nullptr)};
  int startTime {stoi(string(argv[++argi]), nullptr, 10)};
  int endTime {stoi(string(argv[++argi]), nullptr, 10)};
  int timeInc {stoi(string(argv[++argi]), nullptr, 10)};
  string posFile (argv[++argi]);
  string outFile (argv[++argi]);
  
  PositionReader reader;
  reader.open(posFile, numOfBeads, lx, ly, lz, timeInc);
  if (!reader.isOpen()){
    cout << "Problem with reading the position file!" << endl;
    return 1;
  }
  
  ofstream writer;
  writer.open(outFile);
  if (!writer){
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }
  writer << std::setprecision(5) << std::fixed;  
  
  int t, time {};
  double z, wallPos {lz/2.0};
  vector<int> wallBeadCount (typesOfBeads, 0);
  vector<int> beadCount (typesOfBeads, 0);
  while (reader.nextFrame()){
    time = reader.getTime();
    if (time >= startTime && time <= endTime){
      std::fill(wallBeadCount.begin(),wallBeadCount.end(),0);
      std::fill(beadCount.begin(),beadCount.end(),0);
      for (int i {}; i < numOfBeads; i++){
	z = reader.getUnwrappedPosition(i, 2);
	t = reader.getType(i);
	if (z > (wallPos-wallDist)){
	  wallBeadCount[t-1]++;
	}
	beadCount[t-1]++;
      }
      writer << time;
      for (int i {}; i < typesOfBeads; i++){
	writer << " " << 
	  static_cast<double>(wallBeadCount[i])/
	  static_cast<double>(beadCount[i]);
      }
      writer << endl;
    }
  }

  reader.close();
  writer.close();
}
