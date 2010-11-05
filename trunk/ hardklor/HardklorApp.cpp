#include <fstream>
#include <string>

#include "CHardklor.h"
#include "CHardklorParser.h"
#include "CHardklorSetting.h"

using namespace std;

int main(int argc, char* argv[]) {
  int i;
  bool bConf;
  char tstr[512]="\0";
  fstream fptr;

	CAveragine *averagine;
	CMercury8 *mercury;

	if(argc==1){
		cout << "Hardklor v1.34, August 18, 2010\n";
		cout << "Usage:\t\thardklor <MS1 file> <output file> [parameters]\n";
		cout << "\t\thardklor -conf <config file>\n";
		cout << "Parameters:\tSee documentation" << endl;
		exit(1);
	}
  
  bConf = false;
  for(i=0;i<argc;i++){
    if(strcmp(argv[i],"-conf")==0) {
      bConf = true;
      break;
    };
  };
  
  
  CHardklorParser hp;
  if(bConf) {
    hp.parseConfig(argv[i+1]);
  } else {
    for(i=1;i<argc;i++){
      strcat(tstr," ");
      strcat(tstr,argv[i]);
    };
    hp.parse(tstr);
  };
  

  //Create all the output files that will be used
  for(i=0;i<hp.size();i++){
    fptr.clear();
    fptr.open(&hp.queue(i).outFile[0],ios::out);
    fptr.close();
  };

	averagine = new CAveragine(hp.queue(0).MercuryFile,hp.queue(0).HardklorFile);
	mercury = new CMercury8(hp.queue(0).MercuryFile);

  CHardklor h(averagine,mercury);

  for(i=0;i<hp.size();i++) {
    h.GoHardklor(hp.queue(i));
  };

	delete averagine;
	delete mercury;
  
  return 0;
  
};
