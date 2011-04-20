#include "CPeriodicTable.h"

using namespace std;

/*
CPeriodicTable::CPeriodicTable(){
  table = new vector<element>;
  loadTable("Hardklor.dat");
};
*/

CPeriodicTable::CPeriodicTable(char* c){
  table = new vector<element>;
  loadTable(c);
};

CPeriodicTable::~CPeriodicTable(){
  delete table;
};

element& CPeriodicTable::at(int i){
  return table->at(i);
};

void CPeriodicTable::loadTable(char* c){
  FILE* fptr;
  element e;
  
  fptr = fopen(c,"rt");
  if(fptr==NULL) {
    cout << "Cannot open periodic table!" << endl;
    return;
  };
  
  /* 
     This loop reads in table entries, line by line.
     It has basic error detection (missing data), but additional
     checks should be made to confirm the content of data
  */
  while(!feof(fptr)){

		fscanf(fptr,"%d\t%s\t%lf\n",&e.atomicNum,e.symbol,&e.mass);   
    table->push_back(e);
    
  };

  fclose(fptr);
  
};

int CPeriodicTable::size(){
  return table->size();
};
