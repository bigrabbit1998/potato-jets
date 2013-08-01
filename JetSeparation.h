#ifndef JETMATCHING__HH
#define JETMATCHING__HH 1

#include "TLorentzVector.h"
#include "fjClustering.h"
#include "MyEvent.h"

class JetMatching{

 public:
  
  typedef vector<fastjet::PseudoJet> Pseudojets;
  typedef vector<Pythia8::Particle> Particlejets;
  
  JetMatching();
  
  
  /*Jets(void);
  Jets(const Jets &src);
  void operator=(const Jets &src);
  */

  //jet selection cuts
  Pseudojets TrimJets(double, double, Pseudojets&);

  Pseudojets OverlapRemoval(const  Particlejets&,const Pseudojets&);
  /* Jets OverlapRemoval(const Jets &input_jets, const Pythia8:: Particles &ptcls) {
    Jets myjets=input_jets;
    for (size_t i=0;i<ptcls.size();++i) {
      myjets=OverlapRemoval(myjets,ptcls[i]);
    }
    return myjets;
    }*/
  
  //matching is done here
  void SetParam( double,double,double);
  Pseudojets Match(const Particlejets&,const Pseudojets&);
  bool ComparePt(fastjet::PseudoJet, fastjet::PseudoJet );
  bool ChiSquare(TLorentzVector, TLorentzVector );
  void PrintMatches();
  Pseudojets RemoveSubset(const Pseudojets &,const Pseudojets& );


 private:
  
  double m_DeltaR, m_ptmin, m_etamax;
  
  
  
};

#endif
