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
#ifndef _HARDKLORTYPES_H
#define _HARDKLORTYPES_H

#include "MSToolkitTypes.h"

#include <string>
#include <vector>

enum specType {
	OrbiTrap,
	TOF,
	QIT,
	FTICR
};

enum hkAlgorithm {
	Basic,
	SemiComplete,
	SemiCompleteFast,
	Dynamic,
	DynamicSemiComplete,
	SemiSubtractive,
	FewestPeptides,
	FewestPeptidesChoice,
	FastFewestPeptides,
	FastFewestPeptidesChoice,
	Version2
};

enum hkOutputFormat{
	hk,
	pepXML
};

typedef struct {
  char atom[3];
  int isotope;
  double ape;
} sEnrich;

typedef struct {
  std::string molecule;
  int iLower;
  int iUpper;
} sMolecule;

typedef struct {
  int iLower;
  int iUpper;
} sInt;

typedef struct{
  double dLower;
  double dUpper;
} sDouble;

typedef struct{
  float fLower;
  float fUpper;
} sFloat;

typedef struct{
  int iValue;
  double dValue;
} sID;

typedef struct{
  int atomNum;
  int isotope;
  double ape;
} sEnrichMercury;

typedef struct{
  double mz;
  float intensity;
  int index;
} sSplit;

//Do this better
enum ScanType{
  Zoom,
  UltraZoom,
  IonSpec2,
  Other
};

typedef struct {
  char id[5];
  double mz;
  double monoMass;
  double shft;
  double abun;
  int charge;
  char seq[31];
  int C;
  int H;
  int O;
  int N;
  int S;
  std::vector<sID> *enrich;
} peps;



//for Hardklor diagnostics
typedef struct hkDPeakMatch {
  double modelPeak = 0;
  double obsPeak = 0;
  double ppm = -100; //not necessary to store this, but feeling lazy
} hkDPeakMatch;
typedef struct hkDComparison {
  int count = 0;  //number of peaks in use
  int match = 0;  //matches in the final result
  int maxMatch = 0;
  double firstCorr = 0;
  double lastCorr = 0;
  hkDPeakMatch peaks[30];  //should be more than necessary. ever.
} hkDComparison;




typedef struct pepHit{
	int basePeakIndex;
	int charge;
	int lowIndex;
	int highIndex;
	int variantIndex;
  int matchCount;          //new
  int consMatchCount;      //new
  int modelSize;           //new
	float intensity;
	float area;
	double massShift;
	double lowMZ;
	double highMZ;
	double monoMass;
	double corr;
  double modelArea;        //new
	char averagine[32];
  hkDComparison diagnostic; //new
} pepHit;

typedef struct mercuryModel{
	float area;
	int size;
	double zeroMass;
	MSToolkit::Peak_T* peaks;
} mercuryModel;

//for storing results to memory for modular Hardklor runs
typedef struct hkMem{
  int charge;
  int scan;
  float intensity; //or area depending on settings.
  double corr;
  double monoMass;
  double mz;
  char mods[32];
} hkMem;

#endif
