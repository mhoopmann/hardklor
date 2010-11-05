#include "CNoiseReduction.h"

CNoiseReduction::CNoiseReduction(MSReader* msr, CHardklorSetting& hs){
  r=msr;
  cs=hs;
  s = new Spectrum* [cs.noiseWindow];
  for(int i=0;i<cs.noiseWindow;i++) s[i]=new Spectrum;
  pos=0;
}

CNoiseReduction::~CNoiseReduction(){
  for(int i=0;i<cs.noiseWindow;i++) delete s[i];
  delete s;
  s=NULL;
}

//Calculates the resolution (FWHM) of a peak
double CNoiseReduction::calcFWHM(double mz){
	double deltaM;
	switch(cs.msType){
	case OrbiTrap:
		deltaM = mz * sqrt(mz) / (20*cs.res400); //(20*userParams.res400);  //sqare root of 400
		break;
	case TOF:
		deltaM = mz / cs.res400;
		break;
	case QIT:
		deltaM = cs.res400;
		break;
	case FTICR:
	default:
		deltaM = mz * mz / (400*cs.res400);
		break;
	}
	return deltaM;
}

bool CNoiseReduction::DeNoise(Spectrum& sp){

  Spectrum* tmpSpec;
  double ppm;
  int i,j,k;
  int index;
  int center;
  int matchCount;
  char cFilter1[256];
  char cFilter2[256];

  sp.clear();

  //Load data the first time
  if(pos==0){
    if((cs.scan.iLower>0)) {
      r->readFile(&cs.inFile[0],*s[0],cs.scan.iLower);
      for(i=1;i<cs.noiseWindow;i++){
        r->readFile(NULL,*s[i]);
        if(s[i]->getScanNumber()==0){
          cs.noiseWindow=i;
          if(cs.noiseWindow==1){
            cout << "Warning: Persistent Peaks algorithm cannot be used on a single spectrum. Skipping analysis." << endl;
            return false;
          }
          break;
        }
      }
    } else {
      r->readFile(&cs.inFile[0],*s[0]);
      for(i=1;i<cs.noiseWindow;i++){
        r->readFile(NULL,*s[i]);
        if(s[i]->getScanNumber()==0){
          cs.noiseWindow=i;
          if(cs.noiseWindow==1){
            cout << "Warning: Persistent Peaks algorithm cannot be used on a single spectrum. Skipping analysis." << endl;
            return false;
          }
          break;
        }
      }
    }

    //Find peaks (if not reading centroid data)
    if(!cs.centroid){
      for(i=0;i<cs.noiseWindow;i++){
        //Assume High resolution data at all times
        FirstDerivativePeaks(*s[i],1);
      }
    }
  }

  //Initialize data
  center=cs.noiseWindow/2;

  for(j=pos;j<cs.noiseWindow/2;j++){
    pos++;
    for(i=0;i<s[j]->size();i++){
      //Find closest peak
      matchCount=1;
      for(k=0;k<cs.noiseWindow;k++){
        if(k==j) continue;
        if(s[k]->size()<1) continue;
        s[j]->getRawFilter(cFilter1,256);
        s[k]->getRawFilter(cFilter2,256);
        if(strcmp(cFilter1,cFilter2)!=0) continue;
        index = NearestPeak(*s[k],s[j]->at(i).mz);
        ppm=fabs( (s[k]->at(index).mz-s[j]->at(i).mz)/s[j]->at(i).mz * 1000000);
        if(ppm<cs.ppm) matchCount++;
      }
      if(matchCount>=cs.noiseMatch) sp.add(s[j]->at(i));

    }
    sp.setScanNumber(s[j]->getScanNumber());
    sp.setScanNumber(s[j]->getScanNumber(true),true);
    sp.setRTime(s[j]->getRTime());
    return true;
  }

  tmpSpec = new Spectrum;
  r->readFile(NULL,*tmpSpec);
  while(tmpSpec->getScanNumber()>0){

    for(i=0;i<s[center]->size();i++){

      //Find closest peak
      matchCount=1;
      for(k=0;k<cs.noiseWindow;k++){
        if(k==center) continue;
        if(s[k]->size()<1) continue;
        s[center]->getRawFilter(cFilter1,256);
        s[k]->getRawFilter(cFilter2,256);
        if(strcmp(cFilter1,cFilter2)!=0) continue;
        index = NearestPeak(*s[k],s[center]->at(i).mz);
        ppm=fabs( (s[k]->at(index).mz-s[center]->at(i).mz)/s[center]->at(i).mz * 1000000);
        if(ppm<cs.ppm) matchCount++;
      }
      if(matchCount>=cs.noiseMatch) sp.add(s[center]->at(i));
        
    }

    sp.setScanNumber(s[center]->getScanNumber());
    sp.setScanNumber(s[center]->getScanNumber(true),true);
    sp.setRTime(s[center]->getRTime());
    
    //shift spectra
    delete s[0];
    for(i=0;i<cs.noiseWindow-1;i++) s[i]=s[i+1];
    s[cs.noiseWindow-1]=tmpSpec;
    if(!cs.centroid) FirstDerivativePeaks(*s[cs.noiseWindow-1],1);

    return true;

  }
  
  for(j=pos;j<cs.noiseWindow;j++){
    pos++;
    for(i=0;i<s[j]->size();i++){
      //Find closest peak
      matchCount=1;
      for(k=0;k<cs.noiseWindow;k++){
        if(k==j) continue;
        if(s[k]->size()<1) continue;
        s[j]->getRawFilter(cFilter1,256);
        s[k]->getRawFilter(cFilter2,256);
        if(strcmp(cFilter1,cFilter2)!=0) continue;
        index = NearestPeak(*s[k],s[j]->at(i).mz);
        ppm=fabs( (s[k]->at(index).mz-s[j]->at(i).mz)/s[j]->at(i).mz * 1000000);
        if(ppm<cs.ppm) matchCount++;
      }
      if(matchCount>=cs.noiseMatch) sp.add(s[j]->at(i));

    }
    sp.setScanNumber(s[j]->getScanNumber());
    sp.setScanNumber(s[j]->getScanNumber(true),true);
    sp.setRTime(s[j]->getRTime());
    return true;
  }

  return false;

}

//First derivative method taken from CSpecAnalyze, returns base peak intensity of the set
void CNoiseReduction::FirstDerivativePeaks(Spectrum& s, int winSize){
  int i,j;
  float maxIntensity;
  int bestPeak;
  bool bLastPos;
  Spectrum gp;

	int nextBest;
	double FWHM;
	Peak_T centroid;

  bLastPos=false;
  for(i=0;i<s.size()-winSize;i++){

    if(s.at(i).intensity<s.at(i+winSize).intensity) {
      bLastPos=true;
      continue;
    } else {
      if(bLastPos){
				bLastPos=false;
	
				//find max and add peak
				maxIntensity=0;
				for(j=i;j<i+winSize;j++){
				  if (s.at(j).intensity>maxIntensity){
				    maxIntensity=s.at(j).intensity;
				    bestPeak = j;
				  }
				}

				//Best estimate of Gaussian centroid
				//Get 2nd highest point of peak
				if(bestPeak==s.size()-1){
					nextBest=bestPeak-1;
				} else if(s.at(bestPeak-1).intensity > s.at(bestPeak+1).intensity){
					nextBest=bestPeak-1;
				} else {
					nextBest=bestPeak+1;
				}

				//Get FWHM
				FWHM = calcFWHM(s.at(bestPeak).mz);

				//Calc centroid MZ (in three lines for easy reading)
				centroid.mz = pow(FWHM,2)*log(s.at(bestPeak).intensity/s.at(nextBest).intensity);
				centroid.mz /= GC*(s.at(bestPeak).mz-s.at(nextBest).mz);
				centroid.mz += (s.at(bestPeak).mz+s.at(nextBest).mz)/2;

				//Calc centroid intensity
				centroid.intensity=s.at(bestPeak).intensity/exp(-pow((s.at(bestPeak).mz-centroid.mz)/FWHM,2)*GC);

				//some peaks are funny shaped and have bad gaussian fit.
				//if error is more than 10%, keep existing intensity
				if( fabs((s.at(bestPeak).intensity - centroid.intensity) / centroid.intensity * 100) > 10 ||
            //not a good check for infinity
            centroid.intensity>999999999999.9 ||
            centroid.intensity < 0 ) {
					centroid.intensity=s.at(bestPeak).intensity;
				}

				//Hack until I put in mass ranges
				if(centroid.mz<0 || centroid.mz>2000) {
					//do nothing if invalid mz
				} else {
					gp.add(centroid);
				}
				i+=winSize-1;
      }

    }
  }
  
  int scanNumber=s.getScanNumber();
  int scanNumber2=s.getScanNumber(false);
  float rTime=s.getRTime();
  s = gp;
  s.setRTime(rTime);
  s.setScanNumber(scanNumber);
  s.setScanNumber(scanNumber2,true);

}

//Binary search to quickly find the nearest peak
int CNoiseReduction::NearestPeak(Spectrum &s, double mz){
  int pivot=0;
  int width=0;
  int lastWidth=0;
  int best=0;
  double dif=9999999.9;
  double d;

  pivot=s.size()/2;
  width=(int)(pivot/2.0+0.5);
  while(width!=lastWidth && pivot<s.size() && pivot>-1){
    d=fabs(s.at(pivot).mz-mz);
    if(d<dif){
      dif=d;
      best=pivot;
    }
    if(s.at(pivot).mz==mz){
      return pivot;
    } else if(s.at(pivot).mz > mz){
      pivot-=width;
      lastWidth=width;
      width=(int)(width/2.0+0.5);
    } else {
      pivot+=width;
      lastWidth=width;
      width=(int)(width/2.0+0.5);
    }
  }
  if(pivot<s.size() && pivot>-1){
    d=fabs(s.at(pivot).mz-mz);
    if(d<dif) best=pivot;
  }
  return best;
}

