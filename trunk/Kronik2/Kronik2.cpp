#include "CKronik2.h"

void usage();

int main(int argc, char* argv[]){

	cout << "\nKrönik version 2.10, Feb 29 2012\nApp: kronik\nAuthor: Michael Hoopmann\nCopyright 2007-2011\n" << endl;
	cout << "\nDescription: Finds peptide isotope distributions (PIDs) that" << endl;
	cout << "             persist in at least 3 of 4 consecutive scans from" << endl;
	cout << "             Hardklör results.\n" << endl;

  if(argc<3){
    usage();
		exit(0);
  }

  CKronik2 k;
	int i;
	double mass1=600.0;
	double mass2=8000.0;
	double ppmTol=10.0;
	float contam=5.0;
	int gapTol=1;
	int matchTol=3;

	if(argc>3){
		for(i=1;i<argc-2;i+=2){
			if(argv[i][0]!='-') {
				cout << "Invalid flag!\n" << endl;
				usage();
				exit(0);
			}
			switch(argv[i][1]){
			case 'p':	ppmTol=atof(argv[i+1]);	break;
			case 'c': contam=(float)atof(argv[i+1]);	break;
			case 'd': matchTol=atoi(argv[i+1]); break;
			case 'g': gapTol=atoi(argv[i+1]);	break;
			case 'm': mass2=atof(argv[i+1]); break;
			case 'n': mass1=atof(argv[i+1]); break;
			default:
				cout << "Invalid flag!\n" << endl;
				usage();
				exit(0);
				break;
			}
		}
	}

	cout << "Input file:\t" << argv[argc-2] << endl;
	cout << "Output file:\t" << argv[argc-1] << endl;
	cout << "Parameters:" << endl;
	cout << "  Persistence set to " << matchTol << " with " << gapTol << " gaps allowed." << endl;
	cout << "  Mass Tolerance: " << ppmTol << " ppm." << endl;
	cout << "  Contaminants set at " << contam << " minutes persistent." << endl;
	cout << "  Maxmimum peptide size allowed: " << mass2 << " daltons." << endl;
	cout << "  Minimum peptide size allowed: " << mass1 << " daltons." << endl;
	cout << "\n" << endl;

  k.setGapTol(gapTol);
  k.setMatchTol(matchTol);
  k.setPPMTol(ppmTol);

  k.processHK(argv[argc-2]);
	cout << "Total Persistent Peptide Isotope Distributions (PID): " << k.size() << endl;
	k.removeContaminants(contam);
	cout << "PIDs after removing contaminants: " << k.size() << endl;
	k.removeMass(mass1,mass2);
	cout << "PIDs over mass range: " << k.size() << endl;

	FILE* f;
  f=fopen(argv[argc-1],"wt");

  //Heading line
	fprintf(f,"First Scan\tLast Scan\tNum of Scans\tCharge\tMonoisotopic Mass\tBase Isotope Peak\t");
	fprintf(f,"Best Intensity\tSummed Intensity\tFirst RTime\tLast RTime\tBest RTime\tBest Correlation\tModifications\n");

  for(i=0;i<k.size();i++){
		fprintf(f,"%d\t%d\t%d\t%d\t%lf\t%lf\t%f\t%f\t%f\t%f\t%f\t%lf\t%s\n",
																																		   k[i].lowScan,
																																		   k[i].highScan,
                                                                       k[i].datapoints,
																																		   k[i].charge,
																																		   k[i].monoMass,
																																		   k[i].basePeak,
																																		   k[i].intensity,
																																		   k[i].sumIntensity,
																																		   k[i].firstRTime,
																																		   k[i].lastRTime,
																																		   k[i].rTime,
																																		   k[i].xCorr,
																																		   k[i].mods);
	}
  fclose(f);

  return 0;
}

void usage(){
	cout << "Usage: kronik [flags] [input file] [output file]" << endl;
	cout << "\nInput files are Hardklör generated results files." << endl;
	cout << "Output files are in tab-delimited ASCII text format." << endl;
	cout << "\nAll flags must be followed by a positive integer." << endl;
	cout << "Example: kronik -p 5 inputfile.txt outputfile.txt" << endl;
	cout << "\nFlags:" << endl;
	cout << "  -c\tSets the threshold for contaminants. If a peptide isotope\n"
		   << "    \tdistribution persists greater than this number of minutes, it\n"
			 << "    \tis excluded from further analysis. Default: 5.0\n" << endl;
	cout << "  -d\tSets the match tolerance. This number specifies the minimum\n"
		   << "    \tconsecutive scans (allowing for gaps) that a peptide is\n"
			 << "    \tobserved across to be considered persistent. Default: 3\n" << endl;
	cout << "  -g\tSets the gap tolerance. Sets the number of scans a peptide\n"
		   << "    \tcan skip and still be considered if seen again in the folloing\n"
			 << "    \tscan. Default: 1\n" << endl;
	cout << "  -m\tSets the maxmimum mass allowed in daltons. If a peptide isotope\n"
		   << "    \tdistribution's mass is greater than this number, it is\n"
			 << "    \texcluded from further analysis. Default: 8000.0\n" << endl;
	cout << "  -n\tSets the minimum mass allowed in daltons. If a peptide isotope\n"
		   << "    \tdistribution's mass is less than this number, it is excluded\n"
			 << "    \tfrom further analysis. Default: 600.0\n" << endl;
	cout << "  -p\tSets the mass accuracy in parts per million (ppm). Peptides in\n"
		   << "    \tadjacent scans must differ in mass by less than this amount to\n"
			 << "    \tbe called the same. Default: 10.0\n" << endl;
	cout << "\nPlease read the README.txt file for more information on Persistent." << endl;
}