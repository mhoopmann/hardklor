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
#ifndef _CSPLITSPECTRUM_H
#define _CSPLITSPECTRUM_H

#include "CHardklorSetting.h"
#include "CSpecAnalyze.h"
#include "Spectrum.h"
#include <vector>

using namespace std;

class CSplitSpectrum {
 public:
  //Constructors & Destructors
  CSplitSpectrum(Spectrum *spec, CHardklorSetting& sett);
	~CSplitSpectrum();
  
  //Fuctions
  int getNumWindows();
	CSpecAnalyze getWindow(int index);
	void IntersectionAnalysis();
	void MakeAnalysis(double winSize);
	void NoSplitAnalysis();
	void OverlappingAnalysis(double gapSize);
	void SinglePassAnalysis(double gapSize);
	void UnionAnalysis();

	void SetAveragine(CAveragine *a);
	void SetMercury(CMercury8 *m);

	void Centroid(Spectrum& s);

	void NewSNPass(double gapSize);


 private:
  //Data members
  Spectrum *wholeSpec;
	CHardklorSetting userParams;

	CSpecAnalyze goodPeaks;
	CAveragine *averagine;
	CMercury8 *mercury;

	vector<CSpecAnalyze> *setA;
	vector<CSpecAnalyze> *setB;
	vector<CSpecAnalyze> *finalAnalysis;

	vector<float> *s2n;

	vector<int> *aIndex;
	vector<int> *bIndex;

};

#endif

