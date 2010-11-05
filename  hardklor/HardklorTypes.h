#ifndef _HARDKLORTYPES_H
#define _HARDKLORTYPES_H

#include "MSToolkitTypes.h"

#include <string>

using namespace std;

enum specType {
	FTICR,
	OrbiTrap,
	TOF,
	QIT
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
	FastFewestPeptidesChoice
};


typedef struct {
  char atom[3];
  int isotope;
  double ape;
} sEnrich;

typedef struct {
  string molecule;
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

#endif
