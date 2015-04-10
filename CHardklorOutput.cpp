#include "CHardklorOutput.h"

bool CHardklorOutput::openFile(char* fileName, hkOutputFormat format){
	
	if(outFile!=NULL) closeFile();

	outFile=fopen(fileName,"wt");
	if(outFile==NULL){
		cout << "Cannot open output file: " << fileName << endl;
		return false;
	}

	//output header
	fprintf(outFile,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
	fprintf(outFile,"<?xml-stylesheet type=\"text/xsl\" href=\"pepXML_std.xsl\"?>");
	fprintf(outFile,"<ms_analysis date=\"\" summary_xml=\"\" xmlns=\"http://regis-web.systemsbiology.net/pepXML\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://sashimi.sourceforge.net/schema_revision/pepXML/pepXML_v117.xsd\">");
  
	return true;
}

bool CHardklorOutput::exportScan(){

	if(outFile==NULL) {
		cout << "No file open for export." << endl;
		return false;
	}

	return true;
}

void CHardklorOutput::closeFile(){
	if(outFile!=NULL) fclose(outFile);
}