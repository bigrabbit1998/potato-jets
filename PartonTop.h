#ifndef PARTONTOP__HH
#def PARTONTOP__HH 1


#include "fjClustering.h"
#include "MyEvent.h"


class PartonTop::public TNamed {
 public:
  
  typedef Pythia8::Particle pythia;
  typedef vector<Pythia8::Particle> pythiajet;

  //constructor
  PartonTop(Pythia8::Particle top);

  bool IsTop(};
  bool IsAntiTop();
  bool isGluon();
  
  //access paricles

  
 private:
  Pythia8::Particle m_parton;
  

  
  


}
