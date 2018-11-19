// DataManager.h

#ifndef DATA_MANAGER_HPP
#define DATA_MANAGER_HPP

#include <vector>
#include <string>

using std::vector;
using std::string;

class DataManager {

private:
  int haploidNum;
  bool isHaploid;
  int numOfChromo;
  int bpPerBead;

public:

  // Constructor
  DataManager(int haploidNum, bool isHaploid, int bpPerBead);

  // Accessor methods
  void setHaploidNum(int num);
  int getHaploidNum();
  void setBPPerBead(int bpPerBead);
  int getBPPerBead();
  int getNumOfChromo();
  void changeToDuploid();
  void changeToHaploid();

  int getChromoNumber(string ch, const string& prefix);
  bool getFracContent(const string& file, 
		      vector< vector<double> >& fracScore,
		      int headerLines, double thresholdScore, 
		      int chrCol, int startCol, int endCol, 
		      int scoreCol, int totalCol);
};

#endif
