/*
Copyright 2007-2016, Michael R. Hoopmann, Institute for Systems Biology
Michael J. MacCoss, University of Washington

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
#ifndef _CSETTINGS_H
#define _CSETTINGS_H

#include <vector>

using namespace std;

enum ScanType{
  Zoom,
  UltraZoom,
  IonSpec2,
  Other
};

typedef struct {
  int code;
  double value;
} varType;

typedef struct {
  char infile[256];
  char outfile[256];
  ScanType st;
} FileName;

class CSettings {
 public:
  bool QAR;
  bool ROC;
  int chlor;
  int selen;
  int smooth;
  int maxPep;
  int maxCharge;
  int scan;
  int maxDepth;
  double S2N;
  double signal;
  double corrThresh;
  double windowLower;
  double windowUpper;
  vector<FileName*> files;
  vector<varType> variants;
  CSettings();
  CSettings(char*);
  char enrichAtom[3];
  double enrich;
  int enrichTimes;
 protected:
 private:
  void readFile(char*);
};



#endif
