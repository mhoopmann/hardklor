#ifndef _CHARDKLOROUTPUT_H
#define	_CHARDKLOROUTPUT_H

#include "HardklorTypes.h"
#include <iostream>

using namespace std;

class CHardklorOutput {
public:

	bool openFile(char* fileName, hkOutputFormat format);
	bool exportScan();
	bool exportPeptide();
	void closeFile();

private:
	FILE* outFile;
	bool bScanOpen;

};

#endif