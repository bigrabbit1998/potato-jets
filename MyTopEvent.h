#ifndef MYTOPEVENT__HH
#define MYTOPEVENT__HH 1

#include "fjClustering.h"
#include "MyEvent.h"
//#include "PartonTop.h"

//this class will recieve all information of top event 
//and reconstruct the top mass. 
class MyTopEvent/*:: public TNamed*/
{

 public:
  
  typedef vector<fastjet::PseudoJet> mypseudojets;
  typedef vector<Pythia8::Particle> myparticlejets;

  //constructor
  MyTopEvent();


  /* Use parton information to determine decay channel
  bool AllHadronic() { return (m_t.hadronicDecay() && m_tbar.hadronicDecay()); }
  bool DiLeponic() { return (m_t.leptonicDecay()&&m_tbar.leptonicDecay()); }
  bool SemiLeptonic() { return !AllHadroic() && !DiLeptonic(); }
  */

  void ClearJets();
  double Recon_Mass1(const mypseudojets &, const mypseudojets&, Pythia8::Vec4 &, Pythia8::Particle & );
  int BestCombination(const mypseudojets &, Pythia8::Vec4 & );
  double ReturnWmass();
  Pythia8::Vec4 LeptonicW(myparticlejets &, myparticlejets &, myparticlejets &, myparticlejets & );
  Pythia8::Vec4 Summation(myparticlejets &);
  vector<Pythia8::Vec4> Changetype(const mypseudojets & );
 private:
  mypseudojets m_tops, m_bjets, m_lightjets, m_extra, m_tbars;
  double m_mass, m_massW;
  
};

#endif
