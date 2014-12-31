#include <fstream>
#include <string>

#include "CHardklor.h"
#include "CHardklor2.h"
#include "CModelLibrary.h"
#include "CHardklorParser.h"
#include "CHardklorSetting.h"
#include "CHardklorVariant.h"

using namespace std;

int main(int argc, char* argv[]) {
  int i;
	unsigned int j;
  char tstr[512]="\0";
  fstream fptr;

	CAveragine *averagine;
	CMercury8 *mercury;
	CModelLibrary *models;

	cout << "Hardklor v2.17, December 31 2014" << endl;
	cout << "Mike Hoopmann, Mike MacCoss\nCopyright 2007-2014\nUniversity of Washington" << endl;
	if(argc!=2){
		cout << "Usage:\t\thardklor <config file>\n";
		cout << "See documentation for instructions to modify and use config files." << endl;
		exit(1);
	}
  
  CHardklorParser hp;
  hp.parseConfig(argv[1]);  

  //Create all the output files that will be used
  for(i=0;i<hp.size();i++){
    fptr.clear();
    fptr.open(&hp.queue(i).outFile[0],ios::out);
    fptr.close();
  }

	averagine = new CAveragine(hp.queue(0).MercuryFile,hp.queue(0).HardklorFile);
	mercury = new CMercury8(hp.queue(0).MercuryFile);
	models = new CModelLibrary(averagine,mercury);

  CHardklor h(averagine,mercury);
	CHardklor2 h2(averagine,mercury,models);
	vector<CHardklorVariant> pepVariants;
	CHardklorVariant hkv;

  for(i=0;i<hp.size();i++) {
		if(hp.queue(i).algorithm==Version2){
			pepVariants.clear();
			if(!hp.queue(i).noBase) pepVariants.push_back(hkv);
			for(j=0;j<hp.queue(i).variant->size();j++)  pepVariants.push_back(hp.queue(i).variant->at(j));

			models->eraseLibrary();
			models->buildLibrary(hp.queue(i).minCharge,hp.queue(i).maxCharge,pepVariants);
			h2.GoHardklor(hp.queue(i));
		} else {
			h.GoHardklor(hp.queue(i));
		}
  }

	delete models;
	delete averagine;
	delete mercury;
  
  return 0;
  
}
