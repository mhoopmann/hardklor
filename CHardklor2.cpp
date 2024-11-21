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

#include "CHardklor2.h"

using namespace std;
using namespace MSToolkit;

CHardklor2::CHardklor2(CAveragine *a, CMercury8 *m, CModelLibrary *lib){
  averagine=a;
  mercury=m;
	models=lib;
	bEcho=true;
  bMem=false;
	PT=NULL;
}

CHardklor2::~CHardklor2(){
	averagine=NULL;
	mercury=NULL;
	models=NULL;
	if(PT!=NULL) {
		PT=NULL;
	}
}

hkMem& CHardklor2::operator[](const int& index){
  return vResults[index];
}

void CHardklor2::Echo(bool b){
	bEcho=b;
}

//Encodes a numerical data array to a [compressed] base64 string
string CHardklor2::encodeBinary(void* arr, bool bFloat, bool bZlib) {
	int i;

	//zlib if requested - note that it is optional to the function, but this shrinker always uses it
	unsigned char* zCompr = NULL;
	uLong len, zLen;
	if (bZlib) {
		void* arrPtr;
		if (!bFloat) {
			len = (uLong)((vector<double>*)arr)->size() * sizeof(double);
			arrPtr = &((vector<double>*)arr)->at(0);
		} else {
			len = (uLong)((vector<float>*)arr)->size() * sizeof(float);
			arrPtr = &((vector<float>*)arr)->at(0);
		}
		zLen = compressBound(len);
		zCompr = (unsigned char*)calloc((uInt)zLen, 1);
		compress(zCompr, &zLen, (const Bytef*)arrPtr, len);
	}

	//convert to base64
	size_t sz64;
	char* arr64 = NULL;
	if (bZlib) sz64 = zLen;
	else {
		if (bFloat) sz64 = ((vector<float>*)arr)->size() * sizeof(float);
		else sz64 = ((vector<double>*)arr)->size() * sizeof(double);
	}

	i = sz64 % 3;
	if (i > 0)sz64 += (3 - i);
	sz64 = sz64 * 4 / 3;

	arr64 = new char[sz64 + 1];
	i = mzParser::b64_encode(arr64, (char*)zCompr, zLen);
	arr64[sz64] = '\0';

	string st = arr64;
	delete[] arr64;

	return st;
}

int CHardklor2::GoHardklor(CHardklorSetting sett, Spectrum* s){
	
	//Member variables
	MSReader r;
	Spectrum curSpec,c;
	vector<int> v;
	FILE* fout=NULL;
	int TotalScans;
	int manyPep, zeroPep, lowSigPep;
	int iPercent;
	int minutes, seconds;
	int i;
	vector<pepHit> vPeps;

	//initialize variables
	cs=sett;
	loadTime=0;
	analysisTime=0;
	TotalScans=0;
	manyPep=0;
  zeroPep=0;
  lowSigPep=0;
  iPercent=0;
	getTimerFrequency(timerFrequency);

  vResults.clear();

	//Prep output
	if(!cs.exportMzML.empty()){
		mzml=new NeoMzMLParser();
		mzml->mzML.version="1.1.0";
		CnmzCv cv;
		cv.id="MS";
		cv.fullName="Proteomics Standards Initiative Mass Spectrometry Ontology";
		cv.version="4.1.12";
		cv.URI="https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo";
		mzml->mzML.cvList.cv.push_back(cv);
		cv.id = "UO";
		cv.fullName = "Unit Ontology";
		cv.version = "09:04:2014";
		cv.URI = "https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo";
		mzml->mzML.cvList.cv.push_back(cv);
		mzml->mzML.cvList.count=2;
		CnmzSoftware sw;
		sw.id="Hardklor";
		sw.version="2.4.0";
		mzml->mzML.softwareList.software.push_back(sw);
		mzml->mzML.softwareList.count=1;
		CnmzInstrumentConfiguration ic;
		ic.id="IC1";
		mzml->mzML.instrumentConfigurationList.instrumentConfiguration.push_back(ic);
		mzml->mzML.instrumentConfigurationList.count=1;
		CnmzDataProcessing dp;
		dp.id="Hardklor_boxcar_spectral_processing";
		CnmzProcessingMethod pm;
		pm.order=0;
		pm.softwareRef="Hardklor";
		dp.processingMethod.push_back(pm);
		mzml->mzML.dataProcessingList.dataProcessing.push_back(dp);
		mzml->mzML.dataProcessingList.count=1;
		mzml->mzML.run.id=cs.exportMzML;
		mzml->mzML.run.defaultInstrumentConfigurationRef="IC1";
		mzml->mzML.run.spectrumList.defaultDataProcessingRef="Hardklor_boxcar_spectral_processing";
	}

	//For noise reduction
	CNoiseReduction nr(&r,cs);

	//Signature
	//if(bEcho) cout << "\n\nHardklor, v2.06, Mike Hoopmann, Mike MacCoss\nCopyright 2007-2012\nUniversity of Washington\n" << endl;

	//Set the periodic table
	if(PT==NULL) PT=averagine->getPT();

	//Ouput file info to user
	if(bEcho){
		if(s==NULL) cout << "Reading from file: " << cs.inFile << endl;
		if(!bMem) cout << "Writing to file: " << cs.outFile << endl;
	}
  if(cs.fileFormat==dunno) {
    cout << "Unknown file format or bad extension." << endl;
    return -1;
  }

  if (!bMem) fout = fopen(cs.outFile, "wt");

	//read a spectrum
	getExactTime(startTime);

	  //Read in the initial spectrum
  r.setFilter(cs.mzXMLFilter);
  r.setRawFilter(cs.rawFilter);
  if(s!=NULL){
    curSpec=*s;
  } else {
	  if(cs.boxcar==0){
			if (cs.scan.iLower > 0) {
				int scn = cs.scan.iLower;
				r.readFile(&cs.inFile[0], curSpec, scn);
				while (curSpec.getScanNumber() == 0) {
					scn++;
					if (cs.scan.iUpper > 0) {
						if (scn > cs.scan.iUpper) break;
					} else {
						if (scn > r.getLastScan()) break;
					}
					r.readFile(NULL, curSpec, scn);
				}
			} else r.readFile(&cs.inFile[0], curSpec);
	  } else {
		  if(cs.boxcarFilter==0){
        if(!nr.DeNoiseD(curSpec)) curSpec.setScanNumber(0);
			  //if(!nr.DeNoise(curSpec)) curSpec.setScanNumber(0);
			  //do something about this...

		  } else {
        if(!nr.DeNoiseC(curSpec)) curSpec.setScanNumber(0);
		  }
    }
  }

	getExactTime(stopTime);
  tmpTime1=toMicroSec(stopTime);
  tmpTime2=toMicroSec(startTime);
  loadTime+=(tmpTime1-tmpTime2);

	//Check that file was read
  if(curSpec.getScanNumber()==0) {
    if(s!=NULL) {
      cout << "Spectrum is invalid." << endl;
      return -2;
    }
    if(cs.scan.iLower>0) cout << cs.inFile << " is invalid, or requested scan number is of incorrect format." << endl;
    else cout << cs.inFile << " is invalid, or contains no spectrum." << endl;
    return -2;
  }

	//Write scan information to output file.
  if(!bMem){
    if (cs.reducedOutput) WriteScanLine(curSpec, fout, 2);
    else if (cs.xml) WriteScanLine(curSpec, fout, 1);
    else WriteScanLine(curSpec, fout, 0);
  } else {
    currentScanNumber = curSpec.getScanNumber();
  }

	//Output progress indicator
	if(bEcho) cout << iPercent;
  
  //While there is still data to read in the file.
  while(true){

		getExactTime(startTime);
		TotalScans++;
		
		//Smooth if requested
		if(cs.smooth>0) SG_Smooth(curSpec,cs.smooth,4);

		//Centroid if needed; notice that this copy wastes a bit of time.
		//TODO: make this more efficient
		if(cs.boxcar==0 && !cs.centroid) Centroid(curSpec,c);
		else c=curSpec;

		//There is a bug when using noise reduction that results in out of order m/z values
		//TODO: fix noise reduction so sorting isn't needed
		if(c.size()>0) c.sortMZ();

		//Analyze
		QuickHardklor(c,vPeps);

		//export results
		for(i=0;i<(int)vPeps.size();i++){
      if(!bMem){
			  if(cs.reducedOutput) WritePepLine(vPeps[i],c,fout,2);
			  else if(cs.xml) WritePepLine(vPeps[i],c,fout,1);
			  else WritePepLine(vPeps[i],c,fout,0);
      } else {
        ResultToMem(vPeps[i],c);
      }
		}

		//Update progress
		if(bEcho){
			if (r.getPercent() > iPercent){
				if(iPercent<10) cout << "\b";
				else cout << "\b\b";
				cout.flush();
				iPercent=r.getPercent();
				cout << iPercent;
				cout.flush();
			}
		}

		getExactTime(stopTime);
    tmpTime1=toMicroSec(stopTime);
    tmpTime2=toMicroSec(startTime);
    analysisTime+=tmpTime1-tmpTime2;
    
    if(s!=NULL) break;

		//Check if any user limits were made and met
		if( (cs.scan.iUpper == cs.scan.iLower) && (cs.scan.iLower != 0) ){
			break;
		} else if( (cs.scan.iLower < cs.scan.iUpper) && (curSpec.getScanNumber() >= cs.scan.iUpper) ){
			break;
		}

		//export spectrum if requested
		if(!cs.exportMzML.empty()){
			CnmzSpectrum sp;
			sp.id="hardklor scan="+to_string(curSpec.getScanNumber());
			CnmzCvParam c;
			c.cvRef="MS";
			if(curSpec.getMsLevel()==1){
				c.accession="MS:1000579";
				c.name="MS1 spectrum";
			} else {
				c.accession = "MS:1000580";
				c.name = "MSn spectrum";
			}
			sp.cvParam.push_back(c);
			c.accession = "MS:1000511";
			c.name = "ms level";
			c.value=to_string(curSpec.getMsLevel());
			sp.cvParam.push_back(c);
			c.value.clear();
			c.accession = "MS:1000130";
			c.name = "positive scan";
			sp.cvParam.push_back(c);
			c.accession = "MS:1000127";
			c.name = "centroid spectrum";
			sp.cvParam.push_back(c);
			if(curSpec.size()>0){
				c.accession = "MS:1000528";
				c.name = "lowest observed m/z";
				c.value = to_string(curSpec[0].mz);
				c.unitCvRef = "MS";
				c.unitAccession = "MS:1000040";
				c.unitName = "m/z";
				sp.cvParam.push_back(c);
				c.accession = "MS:1000527";
				c.name = "highest observed m/z";
				c.value = to_string(curSpec[curSpec.size()-1].mz);
				sp.cvParam.push_back(c);
			}
			
			vector<float> vAbun; //regardless of how abundance was stored, knock it down to floats
			vector<double> vMz;
			for (int a = 0; a < curSpec.size(); a++) {
				vMz.push_back(curSpec[a].mz);
				vAbun.push_back(curSpec[a].intensity);
			}
			if (vMz.size() > 0) {
				sp.defaultArrayLength = (int)vMz.size(); //update full array size
				sp.binaryDataArrayList.emplace_back();
				sp.binaryDataArrayList[0].binaryDataArray.emplace_back();
				sp.binaryDataArrayList[0].binaryDataArray[0].binary.content = encodeBinary(&vMz, false, true);
				sp.binaryDataArrayList[0].binaryDataArray[0].encodedLength = (int)sp.binaryDataArrayList[0].binaryDataArray[0].binary.content.size();
				CnmzCvParam cv;
				cv.cvRef="MS";
				cv.accession="MS:1000523";
				cv.name="64-bit float";
				sp.binaryDataArrayList[0].binaryDataArray[0].cvParam.push_back(cv);
				cv.accession="MS:1000574";
				cv.name="zlib compression";
				sp.binaryDataArrayList[0].binaryDataArray[0].cvParam.push_back(cv);
				cv.accession = "MS:1000514";
				cv.name = "m/z array";
				cv.unitCvRef="MS";
				cv.unitAccession="MS:1000040";
				cv.unitName="m/z";
				sp.binaryDataArrayList[0].binaryDataArray[0].cvParam.push_back(cv);
			}

			//Process the abundance array if it is going to change
			if (vAbun.size() > 0) {
				//sp.defaultArrayLength = (int)vAbun.size(); //update full array size
				sp.binaryDataArrayList[0].binaryDataArray.emplace_back();
				sp.binaryDataArrayList[0].binaryDataArray[1].binary.content = encodeBinary(&vAbun, true, true);
				sp.binaryDataArrayList[0].binaryDataArray[1].encodedLength = (int)sp.binaryDataArrayList[0].binaryDataArray[1].binary.content.size();
				CnmzCvParam cv; 
				cv.cvRef = "MS";
				cv.accession = "MS:1000521";
				cv.name = "32-bit float";
				sp.binaryDataArrayList[0].binaryDataArray[1].cvParam.push_back(cv);
				cv.accession = "MS:1000574";
				cv.name = "zlib compression";
				sp.binaryDataArrayList[0].binaryDataArray[1].cvParam.push_back(cv);
				cv.accession = "MS:1000515";
				cv.name = "intensity array";
				cv.unitCvRef = "MS";
				cv.unitAccession = "MS:1000131";
				cv.unitName = "number of detector counts";
				sp.binaryDataArrayList[0].binaryDataArray[1].cvParam.push_back(cv);
			}
			mzml->mzML.run.spectrumList.spectrum.push_back(sp);
		}

		//Read next spectrum from file.
		getExactTime(startTime);
		if(cs.boxcar==0) {
			r.readFile(NULL,curSpec);
		} else {
			if(cs.boxcarFilter==0){
				//possible to not filter?
        nr.DeNoiseD(curSpec);
			} else {
			//case 5: nr.DeNoise(curSpec); break; //this is for filtering without boxcar
				nr.DeNoiseC(curSpec);
			}
		}

		getExactTime(stopTime);
		tmpTime1=toMicroSec(stopTime);
		tmpTime2=toMicroSec(startTime);
		loadTime+=(tmpTime1-tmpTime2);

		if(curSpec.getScanNumber()!=0){
			//Write scan information to output file.
			if(cs.reducedOutput){
				WriteScanLine(curSpec,fout,2);
			} else if(cs.xml) {
				fprintf(fout,"</Spectrum>\n");
				WriteScanLine(curSpec,fout,1);
			} else {
				WriteScanLine(curSpec,fout,0);
			}
		} else {
			break;
		}
	}

	if(!bMem) fclose(fout);

	if(!cs.exportMzML.empty()){
		mzml->write(cs.exportMzML.c_str());
	}

	if(bEcho) {
		cout << "\n" << endl;
		cout << "  Total number of scans analyzed: " << TotalScans << endl;

		i=(int)timeToSec(loadTime,timerFrequency);
		minutes = (int)(i/60);
		seconds = i - (60*minutes);
		cout << "\nFile access time: " << minutes << " minutes, " << seconds << " seconds." << endl;
		i=(int)timeToSec(analysisTime,timerFrequency);
		minutes = (int)(i/60);
		seconds = i - (60*minutes);
		cout << "Analysis Time:    " << minutes << " minutes, " << seconds << " seconds." << endl;

		if (minutes==0 && seconds==0){
			cout << "IMPOSSIBLE!!!" << endl;
		} else if(minutes <=2){
			cout << "HOLY FRIJOLE!!" << endl;
		} else if(minutes<=5) {
			cout << "Like lightning!" << endl;
		} else if(minutes<=10){
			cout << "That's pretty damn fast!" << endl;
		} else if(minutes<=20){
			cout << "Monkeys calculate faster than that!" << endl;
		} else if(minutes<=30){
			cout << "You should have taken a lunch break." << endl;
		} else if(minutes<=40){
			cout << "Oi! Too freakin' slow!!" << endl;
		} else {
			cout << "You might be able to eek out some better performance by adjusting your parameters." << endl;
		}
	}

	return 1;

}

int CHardklor2::BinarySearch(Spectrum& s, double mz, bool floor){

	int mid=s.size()/2;
	int upper=s.size();
	int lower=0;

	if(mz>=s[s.size()-1].mz) return s.size()-1;
	if(mz<=s[0].mz) return 0;

	while(s[mid].mz!=mz){
		if(lower>=upper) break;
		if(s[mid].mz>mz){
			upper=mid-1;
			mid=(lower+upper)/2;
		} else {
			lower=mid+1;
			mid=(lower+upper)/2;
		}
	}

	if(floor && s[mid].mz>mz && mid>0) return mid-1;
	else if(!floor && s[mid].mz<mz && mid<s.size()) return mid+1;
	else return mid;

}

//Calculates the resolution (FWHM) of a peak
double CHardklor2::CalcFWHM(double mz,double res,int iType){
	double deltaM;
	switch(iType){
	case 0: //Orbitrap
		deltaM = mz * sqrt(mz) / (20*res);  //sqare root of 400
		break;
	case 1: //TOF
		deltaM = mz / res;
		break;
	case 2: //QIT
		deltaM = res / 5000.0;
		break;
	case 3: //FTICR
	default:
		deltaM = mz * mz / (400*res);
		break;
	}
	return deltaM;
}

//First derivative method, returns base peak intensity of the set
void CHardklor2::Centroid(Spectrum& s, Spectrum& out){
  int i,j;
  float maxIntensity;
  int bestPeak;
  bool bLastPos;

	int nextBest;
	double FWHM;
	Peak_T centroid;

	out.clear();

  //Get boundaries of the spectrum. centroids must be within boundaries.
  double minMZ, maxMZ;
  if(s.size()>0){
    minMZ = s[0].mz;
    maxMZ = s[s.size()-1].mz;
  }

  bLastPos=false;
	for(i=0;i<s.size()-1;i++){

    if(s[i].intensity<s[i+1].intensity) {
      bLastPos=true;
      continue;
    } else {
      if(bLastPos){
				bLastPos=false;

				//find max and add peak
				maxIntensity=0;
				for(j=i;j<i+1;j++){
				  if (s[j].intensity>maxIntensity){
				    maxIntensity=s[j].intensity;
				    bestPeak = j;
				  }
				}

				//Best estimate of Gaussian centroid
				//Get 2nd highest point of peak
				if(bestPeak==s.size()){
					nextBest=bestPeak-1;
				} else if(s[bestPeak-1].intensity > s[bestPeak+1].intensity){
					nextBest=bestPeak-1;
				} else {
					nextBest=bestPeak+1;
				}

				//Get FWHM
				FWHM = CalcFWHM(s[bestPeak].mz,cs.res400,cs.msType);

				//Calc centroid MZ (in three lines for easy reading)
				centroid.mz = (FWHM*FWHM*log(s[bestPeak].intensity/s[nextBest].intensity));
				centroid.mz /= (GAUSSCONST*(s[bestPeak].mz-s[nextBest].mz));
				centroid.mz += ((s[bestPeak].mz+s[nextBest].mz)/2);

				//Calc centroid intensity
				centroid.intensity=(float)(s[bestPeak].intensity/exp(-pow((s[bestPeak].mz-centroid.mz)/FWHM,2)*GAUSSCONST));

				//some peaks are funny shaped and have bad gaussian fit.
				//if error is more than 10%, keep existing intensity
				if( fabs((s[bestPeak].intensity - centroid.intensity) / centroid.intensity * 100) > 10 ||
            //not a good check for infinity
            centroid.intensity>999999999999.9 ||
            centroid.intensity < 0 ) {
					centroid.intensity=s[bestPeak].intensity;
				}
				
				//Centroided peaks must fall within spectrum mass range
				if(centroid.mz<minMZ || centroid.mz>maxMZ) {
					//do nothing if invalid mz, but perhaps find a better way to handle this one day.
				} else {
					out.add(centroid);
				}
			
      }

    }
  }

}

//returns whether or not the peak is still valid. true if peak still exists, false if peak was solved already.
bool CHardklor2::CheckForPeak(vector<Result>& vMR, Spectrum& s, int index){
	double dif=100.0;
	double massDif;
	bool match=false;

	int mid=s.size()/2;
	int upper=s.size();
	int lower=0;

	double FWHM=CalcFWHM(vMR[index].mass,cs.res400,cs.msType);

	if(vMR[index].mass>=s[s.size()-1].mz) {
		mid=s.size()-1;
	} else if(vMR[index].mass<=s[0].mz) {
		mid=0;
	} else {

		while(s[mid].mz!=vMR[index].mass){
			if(lower>=upper) break;
			if(s[mid].mz>vMR[index].mass){
				upper=mid-1;
				mid=(lower+upper)/2;
			} else {
				lower=mid+1;
				mid=(lower+upper)/2;
			}
		}
	}

	dif=fabs(s[mid].mz-vMR[index].mass);
	if(mid>0){
		massDif=fabs(s[mid-1].mz-vMR[index].mass);
		if(massDif<dif) {
			dif=massDif;
			mid--;
		}
	} 
	if(mid<s.size()-1){
		massDif=fabs(s[mid+1].mz-vMR[index].mass);
		if(massDif<dif) {
			dif=massDif;
			mid++;
		}
	}

	if(dif<FWHM){
		if(mask[mid].intensity>0.1) return false;
		else return true;
	}
	return false;

}

int CHardklor2::CompareBPI(const void *p1, const void *p2){
  const pepHit d1 = *(pepHit *)p1;
  const pepHit d2 = *(pepHit *)p2;
	if(d1.basePeakIndex<d2.basePeakIndex) return -1;
	else if(d1.basePeakIndex>d2.basePeakIndex) return 1;
  else return 0;
}

double CHardklor2::LinReg(vector<float>& mer, vector<float>& obs){

  int i,sz;
  double sxx=0,syy=0,sxy=0;

  //Cosine angle correlation
  sxy=0;
  sxx=0;
  syy=0;
	sz=(int)mer.size();
  for(i=0;i<sz;i++){
    sxy += (mer[i]*obs[i]);
    sxx += (mer[i]*mer[i]);
    syy += (obs[i]*obs[i]);
  }

  if(sxx>0 && syy>0 && sxy>0) return sxy/sqrt(sxx*syy);
  else return 0;
    
}

bool CHardklor2::MatchSubSpectrum(Spectrum& s, int peakIndex, pepHit& pep){

	int i,k,n;
	size_t varCount;
	size_t v;
	float max=s[peakIndex].intensity;
	int maxMercuryIndex[3];
	vector<int> charges;
	double dif;
	vector<float> obs;
	vector<float> mer;
	vector<int> vMatchIndex;
	vector<float> vMatchPeak;
	vector<Result> vMR;
	Result r;
	double corr;
	double da;

	//keep track of best hits
	double bestCorr;
	double bestDA;
	int bestCharge;
	double bestMass;
	vector<int> bestMatchIndex;
	vector<float> bestMatchPeak;
	int matchCount;
	int bestMatchCount;
	int thisMaxIndex=0;
	int bestVariant;

	mercuryModel* model=NULL;

	double deltaM = CalcFWHM(s[peakIndex].mz,cs.res400,cs.msType);
	QuickCharge(s,peakIndex,charges);

	bestCorr=0.0;
	bestMatchCount=0;

	//Mark number of variants to analyze
	if(cs.noBase) varCount=cs.variant->size();
	else varCount=cs.variant->size()+1;

	//iterate through all charge states
	for(i=0;i<(int)charges.size();i++){

		for(v=0;v<varCount;v++){

			//get model from library
			dif=0;
			model=models->getModel(charges[i],(int)v,s[peakIndex].mz);
			for(k=0; k<model->size; k++) {
				if(model->peaks[k].intensity>dif){
					dif = model->peaks[k].intensity;
					maxMercuryIndex[0]=k;
				}
			}
			if(k==0) maxMercuryIndex[1]=-1;
			else maxMercuryIndex[1]=maxMercuryIndex[0]-1;		//allow right shift
			maxMercuryIndex[2]=maxMercuryIndex[0]+1;				//allow left shift

			//Apply shift and find mz boundaries
			n=0;
			while(n<3){
				
				if(maxMercuryIndex[n]<0) {
					n++;
					continue;
				}

				//Align the mercury distribution (MD) to the observed peak. The MD can shift in
				//either direction for one peak to adjust for system noise. The width of the MD
				//determines the boundaries for correlation to the observed data.
				double lower=5000.0;
				double upper=0.0;
				double shft = s[peakIndex].mz - model->peaks[maxMercuryIndex[n]].mz;
				
				vMR.clear();
				for(k=0; k<model->size; k++) {				
					r.data=model->peaks[k].intensity;
					r.mass=model->peaks[k].mz+shft;
					vMR.push_back(r);
					if(model->peaks[k].intensity>99.999) thisMaxIndex=(int)vMR.size()-1;
					
					if(r.mass<lower) lower=r.mass;
					if(r.mass>upper) upper=r.mass;
				}
				da=model->area;

				//Add a little buffer to the bounds
				lower-=0.1;
				upper+=0.1;

				//Match predictions to the observed peaks and record them in the proper array.
				corr=PeakMatcherB(vMR,s,lower,upper,deltaM/2,peakIndex,matchCount,vMatchIndex,vMatchPeak);
				//cout << "\tMSS: " << s[peakIndex].mz << " " << s[peakIndex].intensity << "\t" << charges[i] << "\t" << matchCount << "\t" << corr << endl;

				if(corr>bestCorr || (corr>cs.corr && corr+0.025*(matchCount-bestMatchCount)>bestCorr) ){
					bestMatchIndex=vMatchIndex;
					bestMatchPeak=vMatchPeak;
					bestMatchCount=matchCount;
					bestCorr=corr;
					bestMass=model->zeroMass+shft*charges[i];
					bestCharge=charges[i];
					bestDA=da;
					bestVariant=(int)v;
				}

				n++;

			}//while

		}//for v (variant)

	}//for i (charge)
	model=NULL;

	//if above threshold, erase peaks.
	if(bestCorr>cs.corr){

		pep.area=(float)bestDA;
		strcpy(pep.averagine,"");
		pep.basePeakIndex=0;
		pep.charge=bestCharge;
		pep.corr=bestCorr;
		pep.highMZ=0;
		pep.lowMZ=0;
		pep.massShift=0;
		pep.monoMass=bestMass;
		pep.intensity=s[peakIndex].intensity;
		pep.variantIndex=bestVariant;

		//mark which peaks contributed to this analysis
		for(k=0;k<(int)bestMatchIndex.size();k++){
			if(bestMatchPeak[k]*max/s[bestMatchIndex[k]].intensity>0.5) s[bestMatchIndex[k]].intensity=-s[bestMatchIndex[k]].intensity;
			else s[bestMatchIndex[k]].intensity-=bestMatchPeak[k]*max;
		}

		return true;
	}

	return false;

}

double CHardklor2::PeakMatcher(vector<Result>& vMR, Spectrum& s, double lower, double upper, double deltaM, int matchIndex, int& matchCount, int& indexOverlap, vector<int>& vMatchIndex, vector<float>& vMatchIntensity){
	//cout << "PeakMatcher " << lower << "-" << upper << endl;
	vMatchIndex.clear();
	vMatchIntensity.clear();

	vector<float> obs;
	vector<float> mer;
				
	bool match;
	bool bMax=false;
	double corr=0.0;
	double dif;
	double massDif;

	matchCount=0;
	indexOverlap=-1;

	int j,k;
	for(k=0;k<(int)vMR.size();k++) {
		if(vMR[k].data>99.9) bMax=true;

		match=false;
		dif=deltaM;
						
		//look left
		j=matchIndex;
		while(j>-1 && s[j].mz>=lower){
			massDif=s[j].mz-vMR[k].mass;
			if(massDif<-deltaM) break;
			if(fabs(massDif)<dif){
				dif=fabs(massDif);
				match=true;
				matchIndex=j;
			}
			j--;
		}

		//look right
		j=matchIndex+1;
		while(j<s.size() && s[j].mz<=upper){
			massDif=s[j].mz-vMR[k].mass;
			if(massDif>deltaM) break;
			if(fabs(massDif)<dif){
				dif=fabs(massDif);
				match=true;
				matchIndex=j;
			}
			j++;
		}
	
		if(!match) {
      //if expected peak is significant (above 50 rel abun) and has no match, match it to 0.
      if(vMR[k].data>50.0) {
        //cout << " xM: " << vMR[k].mass << "\t0" << endl;
        mer.push_back((float)vMR[k].data);
        obs.push_back(0.0f);
        if(bMax) break;
      }
      
		} else {
			mer.push_back((float)vMR[k].data);
      //cout << " xM: " << vMR[k].mass << "\t" << s[matchIndex].mz << "," << s[matchIndex].intensity << "\t" << (vMR[k].mass- s[matchIndex].mz)/ s[matchIndex].mz*1e6; 
			//if((vMR[k].mass - s[matchIndex].mz) / s[matchIndex].mz * 1e6 >10) cout << "!!!";
			//cout << endl;
			if(mask[matchIndex].intensity>0.1 && vMR[k].data>50) {
				if(indexOverlap<0) {
					indexOverlap=matchIndex;
					//cout << " indexOverlap=" << matchIndex<<endl;
				}
			}
			if(s[matchIndex].intensity<0.1) {
				obs.push_back(0.0f);
			} else {
				matchCount++;
				obs.push_back(s[matchIndex].intensity);
			}
			vMatchIndex.push_back(matchIndex);
			vMatchIntensity.push_back((float)vMR[k].data/100.0f);
		}
	}

	if(matchCount<2) corr=0.0;
	else corr=LinReg(mer,obs);

	//for(j=0;j<mer.size();j++){
  //  cout << " M:" << mer[j] << "\t" << "O:" << obs[j] << endl;
	//}
	//cout << " Corr: " << corr << "(" << matchCount << ")" << endl;

  //remove last matched peaks (possibly overlap with other peaks) but only if they are of low abundance.
	int tmpCount=matchCount;
  while(corr<0.90 && matchCount>2 && mer[mer.size()-1]<50.0){
		mer.pop_back();
		obs.pop_back();
		matchCount--;
		double corr2=LinReg(mer,obs);
		//cout << " Old corr: " << corr << "(" << matchCount+1 << ")" << " New corr: " << corr2 << endl;
		if(corr2>corr) {
			corr=corr2;
			tmpCount=matchCount;
		}
	}
	matchCount=tmpCount;
	//cout << "Done PeakMatcher" << endl;
	return corr;
}

double CHardklor2::PeakMatcherB(vector<Result>& vMR, Spectrum& s, double lower, double upper, double deltaM, int matchIndex, int& matchCount, vector<int>& vMatchIndex, vector<float>& vMatchIntensity){

	vMatchIndex.clear();
	vMatchIntensity.clear();

	vector<float> obs;
	vector<float> mer;
				
	bool match;
	bool bMax=false;
	double corr=0.0;
	double dif;
	double massDif;

	matchCount=0;

	int j,k;
	for(k=0;k<(int)vMR.size();k++) {
		if(vMR[k].data>99.9) bMax=true;
		else bMax=false;

		match=false;
		dif=deltaM;
						
		//look left
		j=matchIndex;
		while(j>-1 && s[j].mz>=lower){
			massDif=s[j].mz-vMR[k].mass;
			if(massDif<-deltaM) break;
			if(fabs(massDif)<dif){
				dif=fabs(massDif);
				match=true;
				matchIndex=j;
			}
			j--;
		}

		//look right
		j=matchIndex+1;
		while(j<s.size() && s[j].mz<=upper){
			massDif=s[j].mz-vMR[k].mass;
			if(massDif>deltaM) break;
			if(fabs(massDif)<dif){
				dif=fabs(massDif);
				match=true;
				matchIndex=j;
			}
			j++;
		}
	
		if(!match) {
			break;
		} else {
			mer.push_back((float)vMR[k].data);
			if(s[matchIndex].intensity<0.1) {
				obs.push_back(0.0f);
			} else {
				matchCount++;
				obs.push_back(s[matchIndex].intensity);
			}
			vMatchIndex.push_back(matchIndex);
			vMatchIntensity.push_back((float)vMR[k].data/100.0f);
		}
	}

	if(matchCount<2) corr=0.0;
	else corr=LinReg(mer,obs);

	int tmpCount=matchCount;
	while(corr<0.90 && matchCount>2){
		mer.pop_back();
		obs.pop_back();
		matchCount--;
		double corr2=LinReg(mer,obs);
		if(corr2>corr) {
			corr=corr2;
			tmpCount=matchCount;
		}
	}
	matchCount=tmpCount;

	return corr;
}

double CHardklor2::PeakMatcher2024_01(vector<Result>& vMR, double areaMR, Spectrum& s, double lower, double upper, double deltaM, int matchIndex, int& matchCount, int& indexOverlap, vector<int>& vMatchIndex, vector<float>& vMatchIntensity, int& cMatchCount, double& areaP, hkDComparison* hkd) {
	//cout << "PeakMatcher 2024.01: " << lower << "-" << upper << endl;
	vMatchIndex.clear();
	vMatchIntensity.clear();

	vector<float> obs;
	vector<float> mer;

	bool match;
	bool bMax = false;
	double corr = 0.0;
	double dif;
	double massDif;

	matchCount = 0;
	indexOverlap = -1;
	if(hkd) hkd->count=0;

	int conMatch=0; //consecutive matches
	int conMax=0;   //max consecutive matches
	int j, k;
	for (k = 0; k < (int)vMR.size(); k++) {
		//cout << "\tVMR:" << vMR.size() << "\t" << vMR[k].mass << "\t" << vMR[k].data << endl;
		if (vMR[k].data > 99.9) bMax = true;
		if(hkd) hkd->peaks[hkd->count].modelPeak=vMR[k].mass;

		match = false;
		dif = 10; //deltaM; //switching to ppm

		//look left
		j = matchIndex;
		while (j > -1 && s[j].mz >= lower) {
			massDif = (s[j].mz - vMR[k].mass)/vMR[k].mass*1e6;
			if (massDif < -10) break;
			if (fabs(massDif) < dif) {
				dif = fabs(massDif);
				match = true;
				matchIndex = j;
			}
			j--;
		}
		//cout << "Left: " << match << "\t" << matchIndex << "\t" << s[matchIndex].mz << "\t" << dif <<endl;

		//look right
		j = matchIndex + 1;
		while (j < s.size() && s[j].mz <= upper) {
			massDif = (s[j].mz - vMR[k].mass)/vMR[k].mass*1e6;
			if (massDif > 10) break;
			if (fabs(massDif) < dif) {
				dif = fabs(massDif);
				match = true;
				matchIndex = j;
			}
			j++;
		}
		//cout << "Right: " << match << "\t" << matchIndex << "\t" << s[matchIndex].mz << "\t" << dif << endl;

		if (!match) {
			//if expected peak is significant (above 50 rel abun) and has no match, match it to 0.
			//if (vMR[k].data > 50.0) {
				if (hkd) {
					hkd->peaks[hkd->count].obsPeak=0; 
					hkd->peaks[hkd->count].ppm=-100;
				}
				//cout << " xM: " << vMR[k].mass << "\t0" << endl;
				mer.push_back((float)vMR[k].data);
				obs.push_back(0.0f);
				if (bMax) break;  //if this was the model's base isotope peak, quit now. it MUST match.
			//}
			conMatch=0;

		} else {
			mer.push_back((float)vMR[k].data);
			//cout << " xM: " << vMR[k].mass << "\t" << s[matchIndex].mz << "," << s[matchIndex].intensity << "\t" << (vMR[k].mass - s[matchIndex].mz) / s[matchIndex].mz * 1e6;
			//if ((vMR[k].mass - s[matchIndex].mz) / s[matchIndex].mz * 1e6 > 10) cout << "!!!";
			//cout << endl;
			if (mask[matchIndex].intensity > 0.1 && vMR[k].data > 50) {
				if (indexOverlap < 0) {
					indexOverlap = matchIndex;
					//cout << " indexOverlap=" << matchIndex << endl;
				}
			}
			if (s[matchIndex].intensity < 0.1) {
				obs.push_back(0.0f);
				if (hkd) {
					hkd->peaks[hkd->count].obsPeak = 0; //a ppm without an obs. peak means this peak had been attenuated.
					hkd->peaks[hkd->count].ppm = (vMR[k].mass - s[matchIndex].mz) / s[matchIndex].mz * 1e6;
				}
			} else {
				if (hkd) {
					hkd->peaks[hkd->count].obsPeak = s[matchIndex].mz;
					hkd->peaks[hkd->count].ppm = (vMR[k].mass - s[matchIndex].mz) / s[matchIndex].mz * 1e6;
				}
				matchCount++;
				conMatch++;
				if (conMatch > conMax) conMax = conMatch;
				obs.push_back(s[matchIndex].intensity);
			}
			if(hkd) hkd->count++;
			vMatchIndex.push_back(matchIndex);
			vMatchIntensity.push_back((float)vMR[k].data / 100.0f);
		}
	}

	if (matchCount < 2) corr = 0.0;
	else corr = LinReg(mer, obs);

	if (hkd) {
		hkd->match=matchCount;
		hkd->firstCorr=hkd->lastCorr=corr;
	}

	//for (j = 0; j < mer.size(); j++) {
	//	cout << " M:" << mer[j] << "\t" << "O:" << obs[j] << endl;
	//}
	//cout << " Corr: " << corr << " Matches: " << matchCount << " Consecutive: " << conMax << endl;

	//remove last matched peaks (possibly overlap with other peaks) but only if they are of low abundance.
	//this bit of refinement is questionable, but exists because dealing with low model peaks and noise is tricky.
	/*
	int tmpCount = matchCount;
	while (corr < 0.90 && matchCount>2 && mer[mer.size() - 1] < 50.0) {
		mer.pop_back();
		obs.pop_back();
		matchCount--;
		double corr2 = LinReg(mer, obs);
		cout << " Old corr: " << corr << "(" << matchCount + 1 << ")" << " New corr: " << corr2 << endl;
		if (corr2 > corr) {
			corr = corr2;
			tmpCount = matchCount;
		}
	}
	matchCount = tmpCount;
	*/
	cMatchCount=conMax;
	float tot=0;
	for(size_t a=0;a<mer.size();a++) if(obs[a]>0) tot+=mer[a]; //tally matched peaks
	//cout << tot << " vs. " << areaMR << " = " << tot/areaMR << endl;
	areaP=tot/areaMR;
	//cout << "Done PeakMatcher" << endl;
	return corr;
}

void CHardklor2::QuickCharge(Spectrum& s, int index, vector<int>& v){

	int i,j;
	double dif;
	double rawCh;
	double rawChR;
	int ch;
	int charge[1000];
  float minIntensity=s[index].intensity/4;

	for(i=cs.minCharge;i<=cs.maxCharge;i++) charge[i]=0;

	//check forward
	for(j=index+1;j<s.size();j++){
		if(s[j].intensity<minIntensity) continue;
			
		dif = s[j].mz - s[index].mz;
		if(dif > 1.1) break;
			
		rawCh=1/dif;
		ch = (int)(rawCh+0.5);
		rawChR=rawCh-(int)rawCh;
		if(rawChR>0.2 && rawChR<0.8) continue;
		if(ch<cs.minCharge || ch>cs.maxCharge) continue;
		charge[ch]=1;
	}
  //if no forward charge, exit now.
  bool bMatch=false;
  for(i=cs.minCharge;i<=cs.maxCharge;i++){
    if(charge[i]>0) {
      bMatch=true;
      break;
    }
  }
  if(!bMatch) {
    v.clear();
    return;
  }

	//check backward
	for(j=index-1;j>=0;j--){
    if (s[j].intensity<minIntensity) continue;
			
		dif = s[index].mz - s[j].mz;
		if(dif > 1.1) break;
			
		rawCh=1/dif;
		ch = (int)(rawCh+0.5);
		rawChR=rawCh-(int)rawCh;
		if(rawChR>0.2 && rawChR<0.8) continue;
		if(ch<cs.minCharge || ch>cs.maxCharge) continue;
		charge[ch]=1;
	}

	v.clear();
	for(i=cs.minCharge;i<=cs.maxCharge;i++){
		if(charge[i]>0) v.push_back(i);
	}

}

void CHardklor2::QuickHardklor(Spectrum& s, vector<pepHit>& vPeps) {
	//iterators
	int i,j,k,n,m,x;
	size_t varCount;
	size_t v;

	//tracking spectrum peak intensities
	float maxHeight=9999999999999.9f;
	float max=0.0f;
	float lowPoint=9999999999999.9f;

	//Mercury storage and variables aligning mercury data (including 1 da shifts)
	mercuryModel* model;
	vector<Result> vMR;
	Result r;
	int maxIndex;
	int thisMaxIndex;
	int maxMercuryIndex[3];
	double da;
	double lower;
	double upper;
	double shft;

	//peak variables
	vector<int> charges;
	double deltaM;
	double dif;
	double corr;
	vector<float> obs;
	vector<float> mer;
	vector<int> vMatchIndex;
	vector<float> vMatchPeak;
	vector<int> vMatchIndex2;
	vector<float> vMatchPeak2;
	int matchCount,matchCount2;
	int indexOverlap;
  double top3[3];

	//refinement variables
	bool keepPH;
	pepHit ph2;
	//pepHit bestKeepPH;
	int lowIndex;
	int highIndex;
	bool corr2;
	double corr3;

	//best hit variables
	double bestCorr;
	double bestLow;
	double bestHigh;
	double bestDA;
	int bestCharge;
	double bestMass;
	vector<int> bestMatchIndex;
	vector<float> bestMatchPeak;
	int bestMatchCount;
	pepHit bestPH;
	bool bestKeepPH;
	int bestOverlap;
	int bestLowIndex;
	int bestHighIndex;
	int bestVariant;
	int bestConsMatchCount = 0;
	double bestCorrModifier = 0;
	double bestModelArea = 0;
	int bestModelSize=0;

	//Results
	pepHit ph;

	//Spectrum variables
	Spectrum origSpec=s;
	Spectrum refSpec=s;
	Spectrum tmpSpec;

	//Diagnostics
	hkDComparison comparison;
	hkDComparison bestComparison;


	//cout << "Next spectrum: " << s.getScanNumber() << endl;

	//create mask
	mask.clear();
	for(i=0;i<s.size();i++) mask.add(s[i].mz,0);

	//find lowest intensity;
	for(i=0;i<s.size();i++){
    //printf("%.6lf\t%.1f\n",s[i].mz, s[i].intensity);
		if(s[i].intensity<1) continue; //zero is not allowed as a low point
		if(s[i].intensity<lowPoint) lowPoint=s[i].intensity;
	}

	//clear results vector
	vPeps.clear();

	//Mark number of variants to analyze
	if(cs.noBase) varCount=cs.variant->size();
	else varCount=cs.variant->size()+1;

	//start the loop through all peaks
	while(true){

		//Find most intense peak. Note that sorting is not possible because
		//peaks change in intensity as they are deconvolved. Also it is advantageous
		//to keep peaks in m/z order
		max=0.0f;
		for(i=0;i<s.size();i++){
			if(s[i].intensity<maxHeight && s[i].intensity>max){
				max=s[i].intensity;
				maxIndex=i;
			}
		}

		//stop searching when we reach lowest original point
		//this prevents overfitting with lots of partial noise peaks
		if(max<lowPoint) break;

		//Get the FWHM estimate for the peak we are at.
		deltaM = CalcFWHM(s[maxIndex].mz,cs.res400,cs.msType);

		//Get the charge states. Note that only remaining peaks are used in the estimate.
		//I'm not sure this is best, but it is simpler and faster
		QuickCharge(s,maxIndex,charges);

		//Reset our correlation and matchcount scores. Then iterate through each charge state and find best
		//match to the peaks.
		bestCorr=0.0;
		bestMatchCount=0;
		for(i=0;i<(int)charges.size();i++){

			//cout << "Next Peak: " << s[maxIndex].mz << " " << s[maxIndex].intensity << "\t" << charges[i] << endl;

			//check all variants
			for(v=0;v<varCount;v++){

        //cout << "Variant: " << v << endl;

				//use model library, align to top 3 peaks
				dif=0;
        top3[0]=top3[1]=top3[2]=0;
        maxMercuryIndex[0]=maxMercuryIndex[1]=maxMercuryIndex[2]=-1;
				model=models->getModel(charges[i],(int)v,s[maxIndex].mz);
				//cout << "Model: " << model->area << endl;
				for(k=0; k<model->size; k++) {
          //cout << "i\t" << model->peaks[k].mz << "\t" << model->peaks[k].intensity << endl;
					//if(model->peaks[k].intensity>dif){
          if(model->peaks[k].intensity>top3[0]){
						//dif = model->peaks[k].intensity;
            top3[2]=top3[1];
            top3[1]=top3[0];
            top3[0]=model->peaks[k].intensity;
            maxMercuryIndex[2]=maxMercuryIndex[1];
            maxMercuryIndex[1]=maxMercuryIndex[0];
						maxMercuryIndex[0]=k;
          } else if(model->peaks[k].intensity>top3[1]) {
            top3[2]=top3[1];
            top3[1]=model->peaks[k].intensity;
            maxMercuryIndex[2]=maxMercuryIndex[1];
            maxMercuryIndex[1]=k;
          } else if(model->peaks[k].intensity>top3[2]) {
            top3[2]=model->peaks[k].intensity;
						maxMercuryIndex[2]=k;
          }
				}
				//if(k==0) maxMercuryIndex[1]=-1;
				//else maxMercuryIndex[1]=maxMercuryIndex[0]-1;		//allow right shift
				//maxMercuryIndex[2]=maxMercuryIndex[0]+1;				//allow left shift

				//Test all three positions for the model. Note that if the first peak is the base peak, then
				//no left shift is tested.
				n=0;
				while(n<3){

					//skip the left shift if already at leftmost peak.
					if(maxMercuryIndex[n]<0) {
						n++;
						continue;
					}

          //cout << "ii\tShift #" << n << endl;

					//Align the mercury distribution (MD) to the observed peak. The MD can shift in
					//either direction for one peak to adjust for system noise. The width of the MD
					//determines the boundaries for correlation to the observed data.
					lower=5000.0;
					upper=0.0;
					shft = s[maxIndex].mz - model->peaks[maxMercuryIndex[n]].mz;
					vMR.clear();
					thisMaxIndex=0;
					da=0;

					//use model library
					for(k=0; k<model->size; k++) {
						
						r.data=model->peaks[k].intensity;
						r.mass=model->peaks[k].mz+shft;
						vMR.push_back(r);
						if(model->peaks[k].intensity>99.999) thisMaxIndex=(int)vMR.size()-1;
					
						if(r.mass<lower) lower=r.mass;
						if(r.mass>upper) upper=r.mass;
					}
					da=model->area;

					//Add a little buffer to the m/z boundaries
					lower-=0.1;
					upper+=0.1;

					//Narrow the search to just the area of the spectrum we need
					lowIndex=BinarySearch(s,lower,true);
					highIndex=BinarySearch(s,upper,false);

					//if max peak shifts to already solved peak, skip
					if(!CheckForPeak(vMR,s,thisMaxIndex)){
						n++;
						continue;
					}

					//Match predictions to the observed peaks and record them in the proper array.
					int cMatchCount=0;
					double areaP=0;
					corr=PeakMatcher2024_01(vMR,da,s,lower,upper,deltaM/2,maxIndex,matchCount,indexOverlap,vMatchIndex,vMatchPeak,cMatchCount,areaP,&comparison);
					//cout << "Top corr: " << corr << endl;
					//cout << "ii.i\t" << s[maxIndex].mz << " " << s[maxIndex].intensity << "\t" << charges[i] << "\t" << matchCount << "\t" << corr << "\t" << indexOverlap << "\t" << maxIndex << "\tn" << n << endl;

					//check any overlap with observed peptides. Overlap indicates deconvolution may be necessary.
					//Deconvolution is at best a rough estimate and is not used if it does not improve the correlation
					//scores.
					keepPH=false;
					//if(indexOverlap>-1 /*&& indexOverlap>maxIndex*/){
					if(indexOverlap==-99) { //this guarantees this block is skipped.

            //cout << "iii\tChecking overlap: " << indexOverlap << "\t" << maxIndex << endl;

						//Find overlapping peptide
						for(m=0;m<(int)vPeps.size();m++){
							if(vPeps[m].basePeakIndex==indexOverlap) break;
						}

						//break out subspectrum; this is done using the original spectrum peak heights, not
						//the current peak heights. The peak heights are then adjusted to account for the currently
						//overlapping peptide model.
						tmpSpec.clear();
						x=0;
						int subIndex=-1;
						for(j=vPeps[m].lowIndex;j<=vPeps[m].highIndex;j++){
							
							while(x<(int)vMatchIndex.size() && j>vMatchIndex[x]) x++;

							//generate temporary subspectrum with peak heights reduced for overlapping model
							if(x<(int)vMatchIndex.size() && j==vMatchIndex[x]){
								tmpSpec.add(origSpec[j].mz,origSpec[j].intensity-vMatchPeak[x]*max);
								x++;
							} else {
								tmpSpec.add(origSpec[j]);
							}

							//get the base peak index of the subspectrum
							if(j==indexOverlap) subIndex=tmpSpec.size()-1;
						}

						//Re-Solve subspectrum and see if it has better correlation
						corr2=MatchSubSpectrum(tmpSpec,subIndex,ph2);
						//cout << "iii.i\tCorr2: " << corr2 << "\t" << ph2.corr << "\t" << origSpec[vPeps[m].basePeakIndex].mz << "\t" << vPeps[m].charge << endl;

						//If correlation is better (or close), go back and try the
						//newly adjusted peaks.
						if(corr2 && ph2.corr+0.025>vPeps[m].corr){
							x=0;

							for(j=lowIndex;j<=highIndex;j++){
								if(x<tmpSpec.size() && s[j].mz==tmpSpec[x].mz){
									refSpec[j].intensity=(origSpec[j].intensity+tmpSpec[x].intensity);
									x++;
								} else {
									refSpec[j].intensity=s[j].intensity;
								}
							}

							//solve merged models
							corr3=PeakMatcher(vMR,refSpec,lower,upper,deltaM/2,maxIndex,matchCount2,indexOverlap,vMatchIndex2,vMatchPeak2);
							//cout << "iii.ii\tCorr3: " << s[maxIndex].mz << " " << s[maxIndex].intensity << "\t" << charges[i] << "\t" << matchCount2 << "\t" << corr3 << "\t" << indexOverlap << endl;

							//keep the new model if it is better than the old one.
							if(corr3>corr) {

								corr=corr3;
								vMatchIndex=vMatchIndex2;
								vMatchPeak=vMatchPeak2;
								matchCount=matchCount2;
								//cout << "Refinement happened somewhere: " << corr3 << endl;

								//refine the overlapping one.
								keepPH=true;
								
							} else {

								//it failed, do nothing and move on.
								keepPH=false;

							}
							
						}

					}//if indexOverlap>-1

					//tCorr applies a modifier for having more/fewer peak matches than a different model.
					//the maximum penalty is 0.025, but bonuses are technically infinite.
					//It is possible to get a slight bump by having more matched peaks than the last best model.
					//However, the actually recorded correlation is the cosine angle without the modifier.
          double tCorr;
          if(bestMatchCount==0) tCorr=0;
          else tCorr=0.025*(matchCount-bestMatchCount)/bestMatchCount;
					//cout << "Old best corr: " << bestCorr << "(" << bestMatchCount << ") This corr: " << corr << "," << corr+tCorr << "(" << matchCount << ")" << endl;
					if(/*corr>bestCorr ||*/ (corr>cs.corr && corr+tCorr>bestCorr) ){
						bestMatchIndex=vMatchIndex;
						bestMatchPeak=vMatchPeak;
						bestMatchCount=matchCount;
						bestConsMatchCount=cMatchCount;
						bestCorrModifier=tCorr;
						bestModelSize=vMR.size();
						bestModelArea=areaP;
						bestCorr=corr;
						bestMass=model->zeroMass+shft*charges[i];
						bestCharge=charges[i];
						bestDA=da;
						bestLow=lower;
						bestHigh=upper;
						bestKeepPH=keepPH;
						bestPH=ph2;
						bestOverlap=m;
						bestLowIndex=lowIndex;
						bestHighIndex=highIndex;
						bestVariant=(int)v;
						bestComparison=comparison;
					}

					n++;
				}//while

			}//for v (variants)

		}//for i (charges)

		//if above threshold, erase peaks.
		if(bestCorr>cs.corr){
			ph.area=(float)bestDA;
			ph.basePeakIndex=maxIndex;
			ph.charge=bestCharge;
			ph.corr=bestCorr;
			ph.highMZ=bestHigh;
			ph.intensity=max;
			ph.lowMZ=bestLow;
			ph.massShift=0.0;
			ph.monoMass=bestMass;
			ph.lowIndex=bestLowIndex;
			ph.highIndex=bestHighIndex;
			ph.variantIndex=bestVariant;
			ph.matchCount=bestMatchCount;
			ph.consMatchCount=bestConsMatchCount;
			ph.modelArea=bestModelArea;
			ph.modelSize=bestModelSize;
			ph.diagnostic=bestComparison;
			if(bestKeepPH){
				//cout << "Changed " << bestOverlap << " which was " << vPeps[bestOverlap].monoMass << "\t" << vPeps[bestOverlap].charge << "\t" << vPeps[bestOverlap].corr << "\t" << vPeps[bestOverlap].variantIndex << endl;
				vPeps[bestOverlap].area=bestPH.area;
				vPeps[bestOverlap].intensity=bestPH.intensity;
				vPeps[bestOverlap].corr=bestPH.corr;
				vPeps[bestOverlap].charge=bestPH.charge;
				vPeps[bestOverlap].monoMass=bestPH.monoMass;
				vPeps[bestOverlap].variantIndex=bestPH.variantIndex;
			}
			vPeps.push_back(ph);
			//cout << vPeps.size()-1 << " New top dog: " << ph.monoMass << "\t" << s[ph.basePeakIndex].mz << "\t" << ph.charge << "\t" << ph.corr << "\t" << ph.variantIndex << "\t" << ph.matchCount << "\t" << ph.consMatchCount << "\t" << ph.modelArea << "\t" << ph.modelSize << endl;
			mask[maxIndex].intensity=100.0f;

			for(k=0;k<(int)bestMatchIndex.size();k++){
				if(bestMatchPeak[k]*max/s[bestMatchIndex[k]].intensity>0.5){
					s[bestMatchIndex[k]].intensity=-s[bestMatchIndex[k]].intensity;
				} else {
					s[bestMatchIndex[k]].intensity-=bestMatchPeak[k]*max;
				}
        //cout << "iv\t" << s[bestMatchIndex[k]].mz << " is now " << s[bestMatchIndex[k]].intensity << endl;
			}
		}

		//set new maximum
		maxHeight=max;

	}

	//Sort results by base peak
  //This sort is expensive. Instead, try sorting before exporting to file. Make sure RefineHits below is not
  //order dependent.
	if(vPeps.size()>0) qsort(&vPeps[0],vPeps.size(),sizeof(pepHit),CompareBPI);

	//cout << "Before RefineHits:" << endl;
	//for(size_t a=0;a<vPeps.size();a++){
	//	cout << vPeps[a].monoMass << "\t" << vPeps[a].charge << "\t" << s[vPeps[a].basePeakIndex].mz << "\t" << vPeps[a].corr << "\t" << vPeps[a].variantIndex << endl;
	//}

	//Refine overfitting based on density
	RefineHits(vPeps,origSpec);

}

//Reduces the number of features (cs.depth) per 1 Da window. This removes a lot
//of false hits resulting from jagged tails on really large peaks. Criteria for
//removal is lowest peak intensity
void CHardklor2::RefineHits(vector<pepHit>& vPeps, Spectrum& s){

	unsigned int i;
	int j;
	double lowp,highp;
	bool bRestart=true;
	vector<pepHit> vTmpHit;
	vector<int> vPepMask;
	list<int> vList;
	list<int>::iterator it;

	//generate an index of the peptides to keep or throw away
	for(i=0;i<vPeps.size();i++) vPepMask.push_back(0);

	//iterate through all hits
	for(i=0;i<vPeps.size();i++){

		//skip anything already marked for removal
		if(vPepMask[i]>0)	continue;

		//put a tolerance around each peak
		lowp=s[vPeps[i].basePeakIndex].mz-0.5;
		highp=s[vPeps[i].basePeakIndex].mz+0.5;

		//put the current hit in the list
		vList.clear();
		vList.push_front(i);

		//find all other hits in the tolerance window
		//look left first
		j=i-1;
		while(j>-1){

			//break out when boundary is reached
			if(s[vPeps[j].basePeakIndex].mz<lowp) break;

			//skip anything marked for removal.
			if(vPepMask[j]>0){
				j--;
				continue;
			}

			//add to list from high to low
			for(it=vList.begin();it!=vList.end();it++){
				if(s[vPeps[j].basePeakIndex].intensity > s[vPeps[*it].basePeakIndex].intensity) break;
			}
			vList.insert(it,j);

			j--;
		}

		//look right
		j=i+1;
		while(j<(int)vPeps.size()){

			//break out when boundary is reached
			if(s[vPeps[j].basePeakIndex].mz>highp) break;

			//skip anything marked for removal.
			if(vPepMask[j]>0){
				j++;
				continue;
			}

			//add to list from high to low
			for(it=vList.begin();it!=vList.end();it++){
				if(s[vPeps[j].basePeakIndex].intensity > s[vPeps[*it].basePeakIndex].intensity) break;
			}
			vList.insert(it,j);

			j++;
		}

		//remove the lowest peptides below threshold (user specified depth)
		if((int)vList.size()>cs.depth){
			it=vList.begin();
			for(j=0;j<cs.depth;j++)	it++;
			for(it=it;it!=vList.end();it++)	vPepMask[*it]=1;
		}

	}

	//copy over the keepers
	for(i=0;i<vPeps.size();i++){
		if(vPepMask[i]) continue;
		vTmpHit.push_back(vPeps[i]);
	}
	vPeps.clear();
	for(i=0;i<vTmpHit.size();i++)vPeps.push_back(vTmpHit[i]);

}

void CHardklor2::ResultToMem(pepHit& ph, Spectrum& s){
  int i,j;
  char mods[32];
  char tmp[16];

  hkm.monoMass = ph.monoMass;
  hkm.charge = ph.charge;
  if(cs.distArea) hkm.intensity = ph.area*ph.intensity;
  else hkm.intensity = ph.intensity;
  hkm.scan = currentScanNumber;
  hkm.mz = s[ph.basePeakIndex].mz;
  hkm.corr = ph.corr;

  //Add mods
  if(!cs.noBase) i=ph.variantIndex-1;
	else i=ph.variantIndex;
	strcpy(mods,"");
  if(i<0) {
	  strcat(mods,"_");
	} else {
    for(j=0;j<cs.variant->at(i).sizeAtom();j++){
		  strcat(mods,PT->at(cs.variant->at(i).atAtom(j).iLower).symbol);
		  sprintf(tmp,"%d",cs.variant->at(i).atAtom(j).iUpper);
      strcat(mods,tmp);
		}
		strcat(mods,"_");
		for(j=0;j<cs.variant->at(i).sizeEnrich();j++){
		  sprintf(tmp,"%.2lf",cs.variant->at(i).atEnrich(j).ape);
      strcat(mods,tmp);
			strcat(mods,PT->at(cs.variant->at(i).atEnrich(j).atomNum).symbol);
      sprintf(tmp,"%d_",cs.variant->at(i).atEnrich(j).isotope);
      strcat(mods,tmp);
		}
	}
  strcpy(hkm.mods,mods);
  vResults.push_back(hkm);
}

void CHardklor2::SetResultsToMemory(bool b){
  bMem=b;
}

int CHardklor2::Size(){
  return (int)vResults.size();
}

void CHardklor2::WritePepLine(pepHit& ph, Spectrum& s, FILE* fptr, int format){
  int i,j;

  if(format==0){
		fprintf(fptr,"P\t%.4lf",ph.monoMass);
		fprintf(fptr,"\t%d",ph.charge);
		if(cs.distArea) fprintf(fptr,"\t%.0f",ph.area*ph.intensity);
		else fprintf(fptr,"\t%.0f",ph.intensity);
		fprintf(fptr,"\t%.4lf",s[ph.basePeakIndex].mz);
		fprintf(fptr,"\t%.4lf-%.4lf",s[ph.lowIndex].mz,s[ph.highIndex].mz);
		fprintf(fptr,"\t0.0000");

		//Add mods
		if(!cs.noBase) i=ph.variantIndex-1;
		else i=ph.variantIndex;
		if(i<0) {
			fprintf(fptr,"\t_");
		} else {
			fprintf(fptr,"\t");
			for(j=0;j<cs.variant->at(i).sizeAtom();j++){
				fprintf(fptr,"%s",PT->at(cs.variant->at(i).atAtom(j).iLower).symbol);
				fprintf(fptr,"%d",cs.variant->at(i).atAtom(j).iUpper);
			}
			fprintf(fptr,"_");
			for(j=0;j<cs.variant->at(i).sizeEnrich();j++){
				fprintf(fptr,"%.2lf",cs.variant->at(i).atEnrich(j).ape);
				fprintf(fptr,"%s",PT->at(cs.variant->at(i).atEnrich(j).atomNum).symbol);
				fprintf(fptr,"%d_",cs.variant->at(i).atEnrich(j).isotope);
			}
		}

		fprintf(fptr,"\t%.4lf",ph.corr);

		fprintf(fptr,"\t%d", ph.modelSize);
		fprintf(fptr,"\t%d", ph.matchCount);
		fprintf(fptr,"\t%d", ph.consMatchCount);
		fprintf(fptr,"\t%.2lf", ph.modelArea);

		for(i=0;i<ph.diagnostic.count;i++){
			fprintf(fptr,"\t%.4lf,%.4lf,%.4lf",ph.diagnostic.peaks[i].modelPeak,ph.diagnostic.peaks[i].obsPeak,ph.diagnostic.peaks[i].ppm);
		}

		fprintf(fptr,"\n");

  } else if(format==1) {
		/*
      fptr << "<Peak Mass=\"" << sa.predPep->at(pepID).GetVariant(varID).GetMonoMass() << "\" ";
      fptr << "ChargeState=\"" << sa.predPep->at(pepID).GetVariant(varID).GetCharge() << "\" ";
      if(cs.distArea) fptr << "Area=\"" << sa.predPep->at(pepID).GetIntensity()*sa.predPep->at(pepID).GetVariant(varID).GetArea() << "\" ";
		  else fptr << "Intensity=\"" << sa.predPep->at(pepID).GetIntensity() << "\" ";
			fptr << "MZ=\"" << sa.predPep->at(pepID).GetMZ() << "\" ";
			fptr << "Window=\"" << sa.peaks.at(0).mz << "-" << sa.peaks.at(sa.peaks.size()-1).mz << "\" ";
			fptr << "SN=\"" << sa.S2NCutoff << "\" ";

      //Add mods
      fptr << "Mod=\"";
			for(j=0;j<sa.predPep->at(pepID).GetVariant(varID).GetHKVariant().sizeAtom();j++){
				fptr << PT->at(sa.predPep->at(pepID).GetVariant(varID).GetHKVariant().atAtom(j).iLower).symbol;
				fptr << sa.predPep->at(pepID).GetVariant(varID).GetHKVariant().atAtom(j).iUpper;
			}
			fptr << "_";
			for(j=0;j<sa.predPep->at(pepID).GetVariant(varID).GetHKVariant().sizeEnrich();j++){
				fptr << setiosflags(ios::fixed) << setprecision(2);
				fptr << sa.predPep->at(pepID).GetVariant(varID).GetHKVariant().atEnrich(j).ape;
				fptr << setiosflags(ios::fixed) << setprecision(4);
				fptr << PT->at(sa.predPep->at(pepID).GetVariant(varID).GetHKVariant().atEnrich(j).atomNum).symbol;
				fptr << sa.predPep->at(pepID).GetVariant(varID).GetHKVariant().atEnrich(j).isotope;
				fptr << "_";
			}
      fptr << "\" ";

			fptr << "Score=\"" << obj.corr << "\"/>" << endl;
			*/

		//reduced output
  } else if(format==2){
		fprintf(fptr,"%.4lf",(ph.monoMass+ph.charge*1.007276466)/ph.charge);
		if(cs.distArea) fprintf(fptr,"\t%.0f",ph.area*ph.intensity);
		else fprintf(fptr,"\t%.0f",ph.intensity);
		fprintf(fptr,"\t%d\n",ph.charge);
	}
}

void CHardklor2::WriteScanLine(Spectrum& s, FILE* fptr, int format){

  if(format==0) {
    fprintf(fptr,"S\t%d\t%.4f\t%s",s.getScanNumber(),s.getRTime(),cs.inFile);

		//For Alex Panchaud, special ZS case
		if(s.getFileType()==ZS || s.getFileType()==UZS){
			if(s.sizeZ()>0){
				for(int i=0;i<s.sizeZ();i++) fprintf(fptr,"\t%d,%.6lf",s.atZ(i).z,s.atZ(i).mh);
			}
		} else {

			//otherwise output precursor info if it exists
			if(s.sizeZ()==1){
				fprintf(fptr,"\t%.4lf\t%d\t%.4lf",s.atZ(0).mh-1.00727649,s.atZ(0).z,s.getMZ());
			} else if(s.sizeZ()>1){
				fprintf(fptr,"\t0.0\t0\t%.4lf",s.getMZ());
			} else {
				fprintf(fptr,"\t0.0\t0\t0.0");
			}
		}
    fprintf(fptr,"\n");

		//For XML output
  } else if(format==1){
    fprintf(fptr,"<Spectrum Scan=\"%d\" ",s.getScanNumber());
		fprintf(fptr,"RetentionTime=\"%.4f\" ",s.getRTime()); 
		fprintf(fptr,"Filename=\"%s\"",cs.inFile);
		if(s.getFileType()==ZS || s.getFileType()==UZS){
			if(s.sizeZ()>0){
				for(int i=0;i<s.sizeZ();i++) fprintf(fptr," PeptideSignal%d=\"%d,%.4lf\"",i,s.atZ(i).z,s.atZ(i).mh);
			}
		} else {
			if(s.sizeZ()==1){
				fprintf(fptr," AccMonoMass=\"%.4lf\" PrecursorCharge=\"%d\" PrecursorMZ=\"%.4lf\"",s.atZ(0).mh-1.00727649,s.atZ(0).z,s.getMZ());
			} else if(s.sizeZ()>1){
				fprintf(fptr," AccMonoMass=\"0.0\" PrecursorCharge=\"0\" PrecursorMZ=\"%.4lf\"",s.getMZ());
			} else {
				fprintf(fptr," AccMonoMass=\"0.0\" PrecursorCharge=\"0\" PrecursorMZ=\"0.0\"");
			}
		}
    fprintf(fptr,">\n");

		//For reduced output
	} else if(format==2) {
		fprintf(fptr, "Scan=%d	RT=%.4f\n", s.getScanNumber(),s.getRTime());
	}
}