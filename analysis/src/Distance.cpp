/* Distance.cpp
 * A program that reads the lammpstrj file and compute
 * the radius of gyration as a function of contour length N
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <vector>
#include "PositionReader.hpp"

using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::vector;

double distanceSq(const double r1[3], const double r2[3]);

int main(int argc, char* argv[]){
  if (argc != 10){
    cout << "Usage: Distance [numOfBeads] [lx] [ly] [lz] "
	 << " [startTime] [endTime] [timeInc] [posFile] [outFile]" << endl;
    return 1;
  }

  int argi {};
  int numOfBeads {stoi(string(argv[++argi]), nullptr, 10)};
  double lx {stod(string(argv[++argi]), nullptr)};
  double ly {stod(string(argv[++argi]), nullptr)};
  double lz {stod(string(argv[++argi]), nullptr)};
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
  
  int time, sep;
  vector<double> endToEndDistAvg (numOfBeads-1, 0.0); // Ignore tether bead
  vector<double> gyrRadiusAvg (numOfBeads-1, 0.0); // Ignore tether bead
  vector<long> count (numOfBeads-1, 0); // Ignore tether bead
  double r1[3], r2[3];

  while (reader.nextFrame()){
    time = reader.getTime();
    if (time >= startTime && time <= endTime){
      for (int i {}; i < numOfBeads-1; i++){
	for (int j {}; j < i; j++){
	  // Calculating (squared) distance between bead i and j
	  for (int k {}; k < 3; k++){
	    r1[k] = reader.getUnwrappedPosition(i,k);
	    r2[k] = reader.getUnwrappedPosition(j,k);
	  }
	  sep = abs(i-j);
	  endToEndDistAvg[sep] += distanceSq(r1,r2);
	  
	  // Calculating radius of gyration from bead i to j
	  double rcm[3] {0.0,0.0,0.0};
	  for (int k {j}; k <= i; k++){
	    for (int l {}; l < 3; l++){
	      rcm[l] += reader.getUnwrappedPosition(k,l);
	    }
	  }
	  for (int l {}; l < 3; l++){
	    rcm[l] /= static_cast<double>(sep+1);
	  }
	  double rg {}, diff {};
	  for (int k {j}; k <= i; k++){
	    for (int l {}; l < 3; l++){
	      diff = reader.getUnwrappedPosition(k,l) - rcm[l];
	      rg += diff*diff;
	    }
	  }
	  rg /= static_cast<double>(sep+1);
	  gyrRadiusAvg[sep] += rg;
	  count[sep]++;
	}
      }
    }
  }
  
  reader.close();

  ofstream writer;
  writer.open(outFile);
  if (!writer){
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }
  writer << std::setprecision(5) << std::fixed;

  for (int i {}; i < numOfBeads-1; i++){
     endToEndDistAvg[i] /= static_cast<double>(count[i]);
     gyrRadiusAvg[i] /= static_cast<double>(count[i]);
     writer << i << " " << endToEndDistAvg[i] 
	    << " " << gyrRadiusAvg[i] << endl;
  }
  
  writer.close();
}

double distanceSq(const double r1[3], const double r2[3]){
  double distSq {}, diff {};
  for (int i {}; i < 3; i++){
    diff = r1[i] - r2[i];
    distSq += diff*diff;
  }
  return distSq;
}
