#include "CHardklorParser.h"

using namespace std;

CHardklorParser::CHardklorParser(){
  vQueue = new vector<CHardklorSetting>;
}

CHardklorParser::~CHardklorParser(){
  delete vQueue;
}

//Takes a command line and breaks it into tokens
//Tokens are then read and used to set global and local parameters
void CHardklorParser::parse(char* cmd) {

  int j;
  bool isGlobal=true;
  sMolecule m;
  char *tok;
  char tmpstr[256];
	char upStr[64];
  string tstr;
  vector<string> vs;

  /*
	//For modifications
	CHardklorVariant v;
	CPeriodicTable* PT;
	int j,k;
	string atom;
	string isotope;
	string percent;
	bool badMod;
	int atomNum;
	bool bNew;
  */

	bool bFile;
	char param[32];

  CHardklorSetting hs;

	//Replace first # with a terminator
	tok=strstr(cmd,"#");
	if(tok!=NULL) strncpy(tok,"\0",1);

	//if we have only white space, exit here
	strcpy(tmpstr,cmd);
	tok=strtok(tmpstr," \t\n");
	if(tok==NULL) return;

	//Check if we have a parameter (has '=' in it) or a file request.
	tok=strstr(cmd,"=");
	if(tok==NULL) bFile=true;
	else bFile=false;

	//Read file and return, if needed
	if(bFile){

		hs=global;
		       
		//on systems that allow a space in the path, require quotes (") to capture
    //complete file name
    strcpy(tmpstr,cmd);

    //Check for quote
    if(tmpstr[0]=='\"'){
			
			//continue reading tokens until another quote is found
			j=1;
			while(true){
				if(j==strlen(tmpstr)){
					cout << "Invalid input file." << endl;
					exit(-1);
				}
				if(tmpstr[j]=='\"') break;
			}
			tmpstr[j]='\0';
			strcpy(hs.inFile,&tmpstr[1]);
			j++;
		} else {
			tok=strtok(tmpstr," \t\n");
			strcpy(hs.inFile,tmpstr);
			j=strlen(tmpstr);
		}

		//Find first non-whitespace
		while(true){
			if(j>=(int)strlen(cmd)){
				cout << "Invalid output file." << endl;
				exit(-1);
			}
			if(cmd[j]!=' ' && cmd[j]!='\t') break;
			j++;
		}

		strcpy(tmpstr,&cmd[j]);

    //Check for quote
    if(tmpstr[0]=='\"'){
			
			//continue reading tokens until another quote is found
			j=1;
			while(true){
				if(j==strlen(tmpstr)){
					cout << "Invalid output file." << endl;
					exit(-1);
				}
				if(tmpstr[j]=='\"') break;
			}
			tmpstr[j]='\0';
			strcpy(hs.outFile,&tmpstr[1]);
			j++;
		} else {
			tok=strtok(tmpstr," \t\n");
			strcpy(hs.outFile,tmpstr);
			j=strlen(tmpstr);
		}

		//cout << hs.inFile << "\t" << hs.outFile << endl;

		hs.fileFormat = getFileFormat(hs.inFile);
		vQueue->push_back(hs);
		return;
	}

	//Read parameter
	tok=strtok(cmd," \t=\n");
	if(tok==NULL) return;
	strcpy(param,tok);
	tok=strtok(NULL," \t=\n");
	if(tok==NULL) {
		warn(param,0);
		return;
	}

	//process parameter
	if(strcmp(param,"algorithm")==0){
		for(j=0;j<(int)strlen(tok);j++) upStr[j]=toupper(tok[j]);
		upStr[j]='\0';
		if(strcmp(upStr,"BASIC")==0) global.algorithm=Basic;
		else if (strcmp(upStr,"VERSION1")==0) global.algorithm=FastFewestPeptides;
		else if (strcmp(upStr,"VERSION2")==0) global.algorithm=Version2;
		else {
			global.algorithm=Version2;
			warn("Unknown algorithm. Defaulting to Version2.",2);
		}

	} else if(strcmp(param,"averagine_mod")==0){
    if(strlen(tok)==1 && tok[0]=='0') global.variant->clear();
    else {
      tstr=tok;
      tok=strtok(NULL," \t\n");
      while(tok!=NULL){
        tstr+=" ";
        tstr+=tok;
        tok=strtok(NULL," \t\n");
      }      
      if(!makeVariant(&tstr[0])) warn("Invalid averagine_mod value. Skipping averagine_mod.",2);
    }

	} else if(strcmp(param,"boxcar_averaging")==0){
		global.boxcar=atoi(tok);
    if(global.boxcar>0 && global.boxcar%2==0) {
			global.boxcar++;
			warn("boxcar_averaging value is even number. Incrementing by 1.",2);
		}

	} else if(strcmp(param,"boxcar_filter")==0){
		global.boxcarFilter=atoi(tok);

	} else if(strcmp(param,"boxcar_filter_ppm")==0){
		global.ppm=atof(tok);

	} else if(strcmp(param,"centroided")==0){
		if(atoi(tok)!=0) global.centroid=true;
		else global.centroid=false;

	} else if(strcmp(param,"charge_algorithm")==0){
		for(j=0;j<(int)strlen(tok);j++) upStr[j]=toupper(tok[j]);
		upStr[j]='\0';
		if(strcmp(upStr,"QUICK")==0) global.chargeMode='Q';
		else if (strcmp(upStr,"FFT")==0) global.chargeMode='F';
		else if (strcmp(upStr,"PATTERSON")==0) global.chargeMode='P';
		else if (strcmp(upStr,"SENKO")==0) global.chargeMode='S';
		else if (strcmp(upStr,"NONE")==0) global.chargeMode='B';
		else {
			global.chargeMode='Q';
			warn("Unknown charge algorithm. Defaulting to Quick.",2);
		}

	} else if(strcmp(param,"charge_max")==0){
		global.maxCharge=atoi(tok);

	} else if(strcmp(param,"charge_min")==0){
		global.minCharge=atoi(tok);

	} else if(strcmp(param,"correlation")==0){
		global.corr=atof(tok);

	} else if(strcmp(param,"depth")==0){
		global.depth=atoi(tok);

	} else if(strcmp(param,"distribution_area")==0){
		if(atoi(tok)!=0) global.distArea=true;
		else global.distArea=false;

	} else if(strcmp(param,"hardklor_data")==0){
		strcpy(global.HardklorFile,tok);

	} else if(strcmp(param,"instrument")==0){
		for(j=0;j<(int)strlen(tok);j++) upStr[j]=toupper(tok[j]);
		upStr[j]='\0';
		if(strcmp(upStr,"FTICR")==0) global.msType=FTICR;
		else if (strcmp(upStr,"ORBITRAP")==0) global.msType=OrbiTrap;
		else if (strcmp(upStr,"TOF")==0) global.msType=TOF;
		else if (strcmp(upStr,"QIT")==0) global.msType=QIT;
		else {
			global.msType=OrbiTrap;
			warn("Unknown instrument type. Defaulting to Orbitrap.",2);
		}

	} else if(strcmp(param,"isotope_data")==0){
    strcpy(global.MercuryFile,tok);

	} else if(strcmp(param,"max_features")==0){
	} else if(strcmp(param,"ms_level")==0){
    if(atoi(tok)==2) global.mzXMLFilter=MS2;
    else global.mzXMLFilter=MS1;

	} else if(strcmp(param,"mz_max")==0){
		global.window.dUpper=atof(tok);

	} else if(strcmp(param,"mz_min")==0){
		global.window.dLower=atof(tok);

	} else if(strcmp(param,"mz_window")==0){
		global.winSize=atof(tok);

	} else if(strcmp(param,"resolution")==0){
		global.res400=atof(tok);

	} else if(strcmp(param,"scan_range_max")==0){
		global.scan.iUpper=atoi(tok);

	} else if(strcmp(param,"scan_range_min")==0){
		global.scan.iLower=atoi(tok);

	} else if(strcmp(param,"sensitivity")==0){
	} else if(strcmp(param,"signal_to_noise")==0){
		global.sn=atof(tok);

	} else if(strcmp(param,"smooth")==0){
		global.smooth=atoi(tok);

	} else if(strcmp(param,"sn_window")==0){
		global.snWindow=atof(tok);

	} else if(strcmp(param,"static_sn")==0){
		if(atoi(tok)!=0) global.staticSN=true;
		else global.staticSN=false;

	} else if(strcmp(param,"xml")==0){
		if(atoi(tok)!=0) global.xml=true;
		else global.xml=false;

	} else {
		warn(param,1);
	}


/*
  //Analyze each token
  for(i=0;i<vs.size();i++){

    //Grab first subtoken, which should be the parameter to change
    tok = strtok(&vs.at(i)[0]," \t\n");
    if(tok==NULL) continue;


		//-ro : reduced output
    } else if(strcmp(tok,"-ro")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) {
				if(isGlobal) global.reducedOutput=true;
			  else hs.reducedOutput=true;
        continue;
      }
			if(strcmp(tok,"true")==0){
				if(isGlobal) global.reducedOutput=true;
			  else hs.reducedOutput=true;
			} else if(strcmp(tok,"false")==0 || atoi(tok)==0) {
				if(isGlobal) global.reducedOutput=false;
			  else hs.reducedOutput=false;
			} else {
				if(isGlobal) global.reducedOutput=true;
			  else hs.reducedOutput=true;
      }


    //-m : Parse averagine molecule variants
    } else if(strcmp(tok,"-m")==0) {
			v.clear();
			PT = new CPeriodicTable(global.HardklorFile);

			while(true){
				tok = strtok(NULL," \t\n");
				if(tok==NULL) break;

				badMod=false;
				if(isdigit(tok[0]) || tok[0]=='.'){
					//we have enrichment
					percent="";
					atom="";
					isotope="";

					//Get the APE
					for(j=0;j<(int)strlen(tok);j++){
						if(isdigit(tok[j]) || tok[j]=='.') {
							if(percent.size()==15){
								cout << "Malformed modification flag: too many digits." << endl;
								badMod=true;
								break;
							}
							percent+=tok[j];
						} else {
							break;
						}
					}
					if(badMod) break;

					//Get the atom
					for(j=j;j<(int)strlen(tok);j++){
						if(isalpha(tok[j])) {
							if(atom.size()==2){
								cout << "Malformed modification flag: invalid atom" << endl;
								badMod=true;
								break;
							}
							atom+=tok[j];
						} else {
							break;
						}
					}
					if(badMod) break;

					//Get the isotope
					for(j=j;j<(int)strlen(tok);j++){
						if(isotope.size()==2){
							cout << "Malformed modification flag: bad isotope" << endl;
							badMod=true;
							break;
						}
						isotope+=tok[j];
					}
					if(badMod) break;

					//format the atom properly
					atom.at(0)=toupper(atom.at(0));
				  if(atom.size()==2) atom.at(1)=tolower(atom.at(1));

					//Get the array number for the atom
					atomNum=-1;
					for(j=0;j<PT->size();j++){
						if(strcmp(PT->at(j).symbol,&atom[0])==0){
							atomNum=j;
							break;
						}
					}

					if(atomNum==-1){
						cout << "Malformed modification flag: Atom not in periodic table" << endl;
						break;
					}

					v.addEnrich(atomNum,atoi(&isotope[0]),atof(&percent[0]));

				} else {
					//we have molecule
					percent="";
					atom="";
					bNew=true;

					//Get the atom
					for(j=0;j<(int)strlen(tok);j++){
      
						//Check for an atom symbol
						if(isalpha(tok[j])) {
	
							//Check if we were working on the count of the previous atom
							if(!bNew) {
								bNew=true;
	  
								//Make sure the atom has uppercase-lowercase letter format;
								atom.at(0)=toupper(atom.at(0));
								if(atom.size()==2) atom.at(1)=tolower(atom.at(1));

								//Look up the new atom
								for(k=0;k<PT->size();k++){
									if(strcmp(PT->at(k).symbol,&atom[0])==0){
										//Add the new atom to the variant
										v.addAtom(k,atoi(&percent[0]));
										break;
									}
								}
	  
								//Clear the fields
								percent="";
								atom="";
							}
	
							//Add this letter to the atom symbol
							if(atom.size()==2){
								cout << "Malformed modification flag: invalid atom" << endl;
								badMod=true;
								break;
							}
							atom+=tok[j];

							if(badMod) break;
	
						} else if(isdigit(tok[j])){
	
							//Whenever we find a digit, we have already found an atom symbol
							bNew=false;
	
							//Add this letter to the atom count
							if(percent.size()==12){
								cout << "Malformed modification flag: unreasonable atom count" << endl;
								badMod=true;
								break;
							}
							percent+=tok[j];

							if(badMod) break;
	
						}      
					}

					//process the last atom
					//Make sure the atom has uppercase-lowercase letter format;
					atom.at(0)=toupper(atom.at(0));
					if(atom.size()==2) atom.at(1)=tolower(atom.at(1));

			    //Look up the new atom
					for(k=0;k<PT->size();k++){
						if(strcmp(PT->at(k).symbol,&atom[0])==0){

							//Add the new atom to the variant
							v.addAtom(k,atoi(&percent[0]));
							break;
	
						}
					}
				}

			}

			if(badMod) continue;
			if(isGlobal) global.variant->push_back(v);
			else hs.variant->push_back(v);

			delete PT;

    //-sl : sensitivity level
    } else if(strcmp(tok,"-sl")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) continue;
      if(atof(tok)>=0) {
				if(isGlobal) global.sl=atoi(tok);
				else hs.sl=atoi(tok);
      }

    //-sna : Signal to Noise Algorithm
    } else if(strcmp(tok,"-sna")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) {
        cout << "Unspecified signal to noise algorithm. Defaulting to STD method." << endl;
        continue;
      }
      if(strcmp(tok,"STD")==0 || strcmp(tok,"std")==0) {
				if(isGlobal) global.sna=0;
				else hs.sna=0;
			} else if (strcmp(tok,"ASN")==0 || strcmp(tok,"asn")==0) {
				if(isGlobal) global.sna=1;
				else hs.sna=1;
      } else if (strcmp(tok,"AVG")==0 || strcmp(tok,"avg")==0) {
				if(isGlobal) global.sna=2;
				else hs.sna=2;
      } else if (strcmp(tok,"MIX")==0 || strcmp(tok,"mix")==0) {
				if(isGlobal) global.sna=3;
				else hs.sna=3;
			} else if (strcmp(tok,"V2")==0 || strcmp(tok,"v2")==0) {
				if(isGlobal) global.sna=4;
				else hs.sna=4;
			} else if (strcmp(tok,"V2ASN")==0 || strcmp(tok,"v2asn")==0) {
				if(isGlobal) global.sna=5;
				else hs.sna=5;
			} else if (strcmp(tok,"V2MIX")==0 || strcmp(tok,"v2mix")==0) {
				if(isGlobal) global.sna=6;
				else hs.sna=6;
			} else {
				cout << "Unkown signal to noise algorithm: " << tok << ". Defaulting to STD method." << endl;
				if(isGlobal) global.sna=0;
				else hs.sna=0;
			}


    //-p : Maximum number of peptide models
    } else if(strcmp(tok,"-p")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) continue;
      if(atoi(tok)>0) {
				if(isGlobal) global.peptide=atoi(tok);
				else hs.peptide=atoi(tok);
      }

    //-mF : mzXML Filter
    } else if(strcmp(tok,"-mF")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) continue;
      if(strcmp(tok,"MS1")==0) {
				if(isGlobal) global.mzXMLFilter=MS1;
				else hs.mzXMLFilter=MS1;
      } else if(strcmp(tok,"MS2")==0){
				if(isGlobal) global.mzXMLFilter=MS2;
				else hs.mzXMLFilter=MS2;
			} else if(strcmp(tok,"MS3")==0){
				if(isGlobal) global.mzXMLFilter=MS3;
				else hs.mzXMLFilter=MS3;
			} else if(strcmp(tok,"ZS")==0){
				if(isGlobal) global.mzXMLFilter=ZS;
				else hs.mzXMLFilter=ZS;
			} else if(strcmp(tok,"UZS")==0){
				if(isGlobal) global.mzXMLFilter=UZS;
				else hs.mzXMLFilter=UZS;
      }
      
    //Reading file or invalid parameter
    } else {
      if(tok[0]=='-') {
	      //we have invalid parameter
        cout << "WARNING: Unknown parameter " << tok << endl;
      } else {
        isGlobal=false;
				hs = global;
        
        //on systems that allow a space in the path, require quotes (") to capture
        //complete file name
        strcpy(tmpstr,tok);

        //Check for quote
        if(tok[0]=='\"'){

          if(tok[strlen(tok)-1]=='\"') {
            strcpy(tmpstr,tok);
          } else {

            //continue reading tokens until another quote is found
            while(true){
              tok=strtok(NULL," \t\n");
              if(tok==NULL) {
                cout << "Invalid input file." << endl;
                exit(-1);
              }
              strcat(tmpstr," ");
              strcat(tmpstr,tok);
              if(tok[strlen(tok)-1]=='\"') break;
            }
          }
          //For some reason, the null string terminator is not places, so place it manually
          strncpy(hs.inFile,&tmpstr[1],strlen(tmpstr)-2);
          hs.inFile[strlen(tmpstr)-2]='\0';

        //If no quote, assume file name has no spaces
        } else {
					strcpy(hs.inFile,tok);
        }

				tok = strtok(NULL," \t\n");
				if(tok==NULL) {
					cout << "Invalid output file." << endl;
					exit(-1);
				}

        //Repeat the entire process for output file
        strcpy(tmpstr,tok);
        if(tok[0]=='\"'){

          if(tok[strlen(tok)-1]=='\"') {
            strcpy(tmpstr,tok);
          } else {

            //continue reading tokens until another quote is found
            while(true){
              tok=strtok(NULL," \t\n");
              if(tok==NULL) {
                cout << "Invalid output file." << endl;
                exit(-1);
              }
              strcat(tmpstr," ");
              strcat(tmpstr,tok);
              if(tok[strlen(tok)-1]=='\"') break;
            }
          }
          //For some reason, the null string terminator is not places, so place it manually
          strncpy(hs.outFile,&tmpstr[1],strlen(tmpstr)-2);
          hs.outFile[strlen(tmpstr)-2]='\0';
          
        } else {
					strcpy(hs.outFile,tok);
        }

				hs.fileFormat = getFileFormat(hs.inFile);
      }
    }
    
  }
	*/

}

//Reads in a config file and passes it to the parser
void CHardklorParser::parseConfig(char* c){
  fstream fptr;
  char tstr[512];

  fptr.open(c,ios::in);
  if(!fptr.good()){
    cout << "Cannot open config file!" << endl;
    return;
  }

  while(!fptr.eof()) {
    fptr.getline(tstr,512);
    if(tstr[0]==0) continue;
    if(tstr[0]=='#') continue;
    parse(tstr);
  }

  fptr.close();
}

CHardklorSetting& CHardklorParser::queue(int i){
  return vQueue->at(i);
}

int CHardklorParser::size(){
  return vQueue->size();
}

//Identifies file format from extension - Must conform to these conventions
MSFileFormat CHardklorParser::getFileFormat(char* c){

	char file[256];
	char ext[256];
	char *tok;

	strcpy(file,c);
	tok=strtok(file,".\n");
	while(tok!=NULL){
		strcpy(ext,tok);
		tok=strtok(NULL,".\n");
	}

	if(strcmp(ext,"ms1")==0 || strcmp(ext,"MS1")==0) return ms1;
	if(strcmp(ext,"ms2")==0 || strcmp(ext,"MS2")==0) return ms2;
	if(strcmp(ext,"bms1")==0 || strcmp(ext,"BMS1")==0) return bms1;
	if(strcmp(ext,"bms2")==0 || strcmp(ext,"BMS2")==0) return bms2;
	if(strcmp(ext,"cms1")==0 || strcmp(ext,"CMS1")==0) return cms1;
	if(strcmp(ext,"cms2")==0 || strcmp(ext,"CMS2")==0) return cms2;
	if(strcmp(ext,"zs")==0 || strcmp(ext,"ZS")==0) return zs;
	if(strcmp(ext,"uzs")==0 || strcmp(ext,"UZS")==0) return uzs;
	if(strcmp(ext,"mzML")==0 || strcmp(ext,"MZML")==0) return mzML;
	if(strcmp(ext,"mzXML")==0 || strcmp(ext,"MZXML")==0) return mzXML;
	if(strcmp(ext,"mgf")==0 || strcmp(ext,"MGF")==0) return mgf;
	if(strcmp(ext,"mz5")==0 || strcmp(ext,"MZ5")==0) return mz5;
  if(strcmp(ext,"raw")==0 || strcmp(ext,"RAW")==0) return raw;
	return dunno;

}

void CHardklorParser::warn(char* c, int i){
	switch(i){
		case 0:
			cout << "Parameter " << c << " has no value." << endl;
			break;
		case 1:
			cout << "Unknown parameter: " << c << endl;
			break;
		case 2:
		default:
			cout << c << endl;
			break;
	}
}

bool CHardklorParser::makeVariant(char* c){

  //For modifications
	CHardklorVariant v;
	CPeriodicTable* PT;
	int j,k;
	string atom;
	string isotope;
	string percent;
	int atomNum;
	bool bNew;
  
  char str[256];
  char* tok;
  strcpy(str,c);

  v.clear();
	PT = new CPeriodicTable(global.HardklorFile);  

  tok=strtok(str," \n");  
  while(tok!=NULL){
    if(isdigit(tok[0]) || tok[0]=='.'){
      //we have enrichment
      percent="";
      atom="";
      isotope="";
      
      //Get the APE
      for(j=0;j<(int)strlen(tok);j++){
        if(isdigit(tok[j]) || tok[j]=='.') {
          if(percent.size()==15){
            warn("Bad averagine_mod: Malformed modification flag, too many digits.",2);
            return false;
          }
          percent+=tok[j];
        } else {
          break;
        }
      }
      
      //Get the atom
      for(j=j;j<(int)strlen(tok);j++){
        if(isalpha(tok[j])) {
          if(atom.size()==2){
            warn("Bad averagine_mod: Malformed modification flag, invalid atom",2);
            return false;
          }
          atom+=tok[j];
        } else {
          break;
        }
      }

      //Get the isotope
      for(j=j;j<(int)strlen(tok);j++){
        if(isotope.size()==2){
          warn("Bad averagine_mod: Malformed modification flag, bad isotope",2);
          return false;
        }
        isotope+=tok[j];
      }

      //format the atom properly
      atom.at(0)=toupper(atom.at(0));
      if(atom.size()==2) atom.at(1)=tolower(atom.at(1));

      //Get the array number for the atom
      atomNum=-1;
      for(j=0;j<PT->size();j++){
        if(strcmp(PT->at(j).symbol,&atom[0])==0){
          atomNum=j;
          break;
        }
      }

      if(atomNum==-1){
        warn("Bad averagine_mod: Malformed modification flag, atom not in periodic table",2);
        return false;
      }

      v.addEnrich(atomNum,atoi(&isotope[0]),atof(&percent[0]));

    } else {
      //we have molecule
      percent="";
      atom="";
      bNew=true;

      //Get the atom
      for(j=0;j<(int)strlen(tok);j++){

        //Check for an atom symbol
        if(isalpha(tok[j])) {

          //Check if we were working on the count of the previous atom
          if(!bNew) {
            bNew=true;

            //Make sure the atom has uppercase-lowercase letter format;
            atom.at(0)=toupper(atom.at(0));
            if(atom.size()==2) atom.at(1)=tolower(atom.at(1));

            //Look up the new atom
            for(k=0;k<PT->size();k++){
              if(strcmp(PT->at(k).symbol,&atom[0])==0){
                //Add the new atom to the variant
                v.addAtom(k,atoi(&percent[0]));
                break;
              }
            }

            //Clear the fields
            percent="";
            atom="";
          }

          //Add this letter to the atom symbol
          if(atom.size()==2){
            warn("Bad averagine_mod: Malformed modification flag, invalid atom",2);
            return false;
          }
          atom+=tok[j];
        
        } else if(isdigit(tok[j])){
	
          //Whenever we find a digit, we have already found an atom symbol
          bNew=false;
          
          //Add this letter to the atom count
          if(percent.size()==12){
            warn("Bad averagine_mod: Malformed modification flag, unreasonable atom count",2);
            return false;
          }
          percent+=tok[j];
        
        }      
      }
      
      //process the last atom
      //Make sure the atom has uppercase-lowercase letter format;
      atom.at(0)=toupper(atom.at(0));
      if(atom.size()==2) atom.at(1)=tolower(atom.at(1));
      
      //Look up the new atom
      for(k=0;k<PT->size();k++){
        if(strcmp(PT->at(k).symbol,&atom[0])==0){
          
          //Add the new atom to the variant
          v.addAtom(k,atoi(&percent[0]));
          break;
        
        }
      }
    }
    tok=strtok(NULL," \n"); 
	}

	global.variant->push_back(v);
	delete PT;
  return true;

}