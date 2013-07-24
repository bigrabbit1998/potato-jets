#ifndef MYHISTO__HH
#define MYHISTO__HH 1

#include <string>
#include<iostream> // needed for io
#include<sstream>  // needed for internal io
#include<vector>
#include <cstdio>
#include "Pythia.h"

class MyHisto {
  
 public:
  MyHisto();
  ~MyHisto();
  
  void Create();
  void Write();
  void CloseFile(); 
  
 private:
  double totalhists;
   
}
#endif





