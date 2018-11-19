/* CombinePosXYZ.cpp
 * A simple code that converts position file to ones compatible
 * with xyz format that is readable by vmd
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::vector;
using std::string;
using std::map;

int main(int argc, char* argv[]){
  // Read input arguments
  if (argc != 5){
    cout << "Usage: CombinePosXYZ [numOfBeads] [numOfOutputBeads] " 
	 << "[posFile] [xyzFile]" << endl;
    return 1;
  }
  
  int numOfBeads {stoi(string(argv[1]), nullptr, 10)};
  int numOfOutputBeads {stoi(string(argv[2]), nullptr, 10)};
  string posFile (argv[3]);
  string xyzFile (argv[4]);

  // xyz element map
  map<int,string> elementMap;
  elementMap[1] = "O";
  elementMap[2] = "N";
  elementMap[3] = "C";
  elementMap[4] = "H";
  elementMap[5] = "F";
  elementMap[6] = "S";

  // For storing bead info
  vector<double> zeroVec(3, 0.0);
  vector< vector<double> >* position 
  {new vector< vector<double> >(numOfBeads, zeroVec)};
  vector<int>* type {new vector<int>(numOfBeads, 1)};
  vector<int> zeroIndex(3, 0);
  vector< vector<int> >* boundaryCount
  {new vector< vector<int> >(numOfBeads, zeroIndex)};

  // Read and convert position file
  string line;
  istringstream iss;
  ifstream reader;
  ofstream writer;
  reader.open(posFile);
  writer.open(xyzFile);

  if (!reader){
    cout << "Problem with opening position file!" << endl;
    return 1;
  }
  if (!writer){
    cout << "Problem wiht opening output file!" << endl;
    return 1;
  }

  double x, y, z;
  int index, t, ix, iy, iz;
  long timestep;

  while (!reader.eof()){
    // Skip header lines
    getline(reader, line);
    getline(reader, line);
    iss.clear();
    iss.str(line);
    iss >> timestep;
    for (int i {}; i < 7 && !reader.eof(); i++){
      getline(reader, line);
    }
    // Read bead position and type
    for (int i {}; i < numOfBeads && !reader.eof(); i++){
      getline(reader, line);
      iss.clear();
      iss.str(line);
      iss >> index >> t >> x >> y >> z >> ix >> iy >> iz;
      (*type)[index-1] = t;
      (*position)[index-1][0] = x;
      (*position)[index-1][1] = y;
      (*position)[index-1][2] = z;
      (*boundaryCount)[index-1][0] = ix;
      (*boundaryCount)[index-1][1] = iy;
      (*boundaryCount)[index-1][2] = iz;
    }
    // Write position
    writer << numOfOutputBeads << endl;
    writer << "Atoms. Timestep: " << timestep << endl;
    for (int i {}; i < numOfOutputBeads; i++){
      writer << elementMap[(*type)[i]] << " "
	     << (*position)[i][0] << " "
	     << (*position)[i][1] << " "
	     << (*position)[i][2] << " "
	     << (*boundaryCount)[i][0] << " "
	     << (*boundaryCount)[i][1] << " "
	     << (*boundaryCount)[i][2] << " "
	     << (*type)[i] << endl;
    }
  }

  reader.close();
  writer.close();

  // Delete resources
  delete position;
  delete type;
  delete boundaryCount;
}
