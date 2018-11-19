// DataManager.cpp

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "DataManager.hpp"

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ifstream;

// Constructor
DataManager::DataManager(int num, bool isHap, int bpb) :
  haploidNum {num}, isHaploid {isHap}, bpPerBead {bpb} {
    if (isHaploid){
      numOfChromo = haploidNum;
    } else {
      numOfChromo = haploidNum*2;
    }
}

// Accessor methods
void DataManager::setHaploidNum(int num){
  haploidNum = num;
  if (isHaploid){
    numOfChromo = haploidNum;
  } else {
    numOfChromo = haploidNum*2;
  }
}

int DataManager::getHaploidNum(){
  return haploidNum;
}

void DataManager::setBPPerBead(int bpb){
  bpPerBead = bpb;
}

int DataManager::getBPPerBead(){
  return bpPerBead;
}

int DataManager::getNumOfChromo(){
  return numOfChromo;
}

void DataManager::changeToDuploid(){
  if (isHaploid){
    isHaploid = false;
    numOfChromo = haploidNum*2;
  }
}

void DataManager::changeToHaploid(){
  if (!isHaploid){
    isHaploid = true;
    numOfChromo = haploidNum;
  }
}

// Convert chromosome key into chromosome number
int DataManager::getChromoNumber(string ch, const string& prefix){
  ch.erase(0, prefix.length());
  if (ch != "" && ch != "Y"){
    if (ch == "X")
      return haploidNum;
    else
      return stoi(ch, nullptr, 10);
  }
  return 0;
}

// Compute the fractional content for each data type
bool DataManager::getFracContent(const string& file, 
				 vector< vector<double> >& fracScore,
				 int headerLines, double thresholdScore, 
				 int chrCol, int startCol, int endCol, 
				 int scoreCol, int totalCol){

  const string chrPrefix {"chr"};
  string token {};
  double resolution = static_cast<double>(bpPerBead);
  double score {}, fracContent {};
  long long start {}, end {}; 
  int chromo {}, startBead {}, endBead {};

  cout << "Reading \"" << file << "\" file ... " << endl;
  
  ifstream reader;
  reader.open(file);
  if (!reader){
	cout << "Unable to read the file. "
	     << "Aborting reading process ..." << endl;
	return false;
  }
  
  // Skip header line
  for (int i {}; i < headerLines; i++)
	getline(reader, token); 
  
  bool reachEOF {false};
  while (!reader.eof()){
    for (int col {}; col < totalCol; col++){
      if (!reader.eof()){
	reader >> token;
	if (col == chrCol)
	  chromo = getChromoNumber(token, chrPrefix);
	else if (col == startCol)
	  start = stol(token, nullptr, 10);
	else if (col == endCol)
	  end = stol(token, nullptr, 10);
	else if (col == scoreCol)
	  score = stod(token, nullptr);
      } else {
	reachEOF = true;
	break;
      }
    }

    if (reachEOF) break;
	
    if (chromo > 0 && score > thresholdScore){
      startBead = start / resolution;
      endBead = end / resolution;
	  
      // For content contained within the same bead
      if (startBead == endBead){
	fracContent = static_cast<double>(end-start)/resolution;
	fracScore[chromo-1][startBead] += fracContent;
	if (!isHaploid){
	  fracScore[chromo+haploidNum-1][startBead] += fracContent;
	}
	// For content spread over multiple beads
      } else {
	// Start bead content
	fracContent = 1.0-(static_cast<double>(start)/resolution
			   - static_cast<double>(startBead));
	fracScore[chromo-1][startBead] += fracContent;
	if (!isHaploid){
	  fracScore[chromo+haploidNum-1][startBead] += fracContent;
	}
	// End bead content
	fracContent = static_cast<double>(end)/resolution
	  - static_cast<double>(endBead);
	fracScore[chromo-1][endBead] += fracContent;
	if (!isHaploid){
	  fracScore[chromo+haploidNum-1][endBead] += fracContent;
	}
	// Other beads in betweeen are completely coded by the content
	for (int i {startBead+1}; i < endBead; i++){
	  fracScore[chromo-1][i] = 1.0;
	  if (!isHaploid){
	    fracScore[chromo+haploidNum-1][i] = 1.0;
	  }
	}
      }
    }
  }
  reader.close();
  return true;
}
