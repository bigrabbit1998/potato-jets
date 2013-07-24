#include "fjClustering.h"
#include "MyEvent.h"
#include "MyTopEvent.h"
//#include "PartonTop.h"

//constructor 
MyTopEvent::MyTopEvent(){
  cout<<"Default constructor called"<<endl;
}

void MyTopEvent::ClearJets() {
  /* m_tops.clear();
  m_tbars.clear();
  m_bjets.clear();
  m_lightjets.clear();
  m_extra.clear();  */
}

/*
bool MYTopEvent:: AllHadronic( ) {

  return (m_t.hadronicDecay() && m_tbar.hadronicDecay());
}

bool MyTopEvent::DiLeponic() { 
  return (m_t.leptonicDecay() && m_tbar.leptonicDecay()); 
}

bool MyTopEvent::SemiLeptonic() {
  return !AllHadroic() && !DiLeptonic(); 
}

*/

double MyTopEvent::Recon_Mass1(const mypseudojets & Bjets, const mypseudojets &Lightjets, Pythia8::Vec4 & rest, Pythia8::Particle& nutrin ){
  double rmass(0),mass(0);
  int bestcombo(0);

  if(Lightjets.size() < 2 || Bjets.size()!=2) return -1.0;

  fastjet::PseudoJet besttwolightjets(0,0,0,0);
  for(size_t k(0); k < 2; ++k){
    besttwolightjets +=(Lightjets[k]);
  }

  vector<Pythia8::Vec4> lightjetsvec4 = Changetype( Lightjets);
  Pythia8::Vec4 hadronicw = operator+(lightjetsvec4[0],lightjetsvec4[1]);

  // m_massW = hadronicw.mCalc();
  m_massW = besttwolightjets.m();

  bestcombo= BestCombination(Bjets, hadronicw);
  vector<fastjet::PseudoJet> remainingjets;
  for(size_t k(0); k < Bjets.size(); ++k)
    if(k != bestcombo)
      remainingjets.push_back(Bjets[k]);

  
  vector<Pythia8::Vec4> lastb= Changetype(remainingjets);

  Pythia8::Vec4 nuuuu(nutrin.px(), nutrin.py(), nutrin.pz(), nutrin.e());
  Pythia8::Vec4 pseudotopm = operator+(lastb[0], rest);
  // pseudotopm=operator+(pseudotopm,nuuuu);

  // remainingjet[0]+=(leptonicw);

  rmass = pseudotopm.mCalc();
  //  rmass=besttwolightjets.m();
  return rmass;
}

vector<Pythia8::Vec4> MyTopEvent::Changetype(const mypseudojets & changetype ){
  vector<Pythia8::Vec4> returnjets;
  for(size_t n(0); n < changetype.size(); ++n){
    double px= changetype[n].px();
    double py= changetype[n].py();
    double pz= changetype[n].pz();
    double E= changetype[n].e();
    Pythia8::Vec4 t(px,py,pz,E);
    returnjets.push_back(t);
  }
  return returnjets;
}


double MyTopEvent::ReturnWmass(){
  return m_massW;
}

//solve for pz after assuming w mass
Pythia8::Vec4 MyTopEvent::LeptonicW(vector<Pythia8::Particle> &nutrinos, vector<Pythia8::Particle> & tau, 
			     vector<Pythia8::Particle> &muon, vector<Pythia8::Particle> &elec ) {
  Pythia8::Vec4 leptonW( 0, 0,0,0); 
  Pythia8::Vec4 nu, tu, mu, el, total;
  tu = Summation(tau);
  mu = Summation(muon);
  el = Summation(elec);
 
  total = operator+(mu,el);
  // total=operator+(total,tau);
  double pwE = 
  
  double pwx = total.px() + nu.px();
  double pwy = total.py() + nu.py();
 
  double pwz = fabs((pwE^2 - pwx^2 - pwy^2 - (80.4)^2)^(1/2));

  
  return leptonW;
}

Pythia8::Vec4 MyTopEvent::Summation(myparticlejets &temporary){
 double px(0),py(0),pz(0),e(0);

  //sum over all leptons
 for(size_t v(0); v < temporary.size(); ++v){
    px+=temporary[v].px();
    py+=temporary[v].py();
    pz+=temporary[v].pz();
    e+= temporary[v].e();
  }
  
  Pythia8::Vec4 sendback(px,py,pz,e);
  return sendback;
}

//returns position of best top-w pair
int MyTopEvent::BestCombination(const mypseudojets & btaggedjets, /* vector<Pythia8::Vec4> &w*/ Pythia8::Vec4 &in ){
  mypseudojets bestbjets;
  
  double tpmass(100);
  double tempdif(0);
  int bpos(0);
  double m_topmass=172.5;

  //choose the first two elements of each vector
  //i.e. the two highest pt jets
  for(size_t a(0) ; a < btaggedjets.size() ; ++a){
    const fastjet::PseudoJet &bf = btaggedjets[a];
    double x=bf.px();
    double y=bf.py();
    double z=bf.pz();
    double E=bf.e();
    Pythia8::Vec4 temc(x,y,z,E);
    
    /* for(size_t v(0); v < w.size(); ++v) {
       const Pythia8::Vec4 &wjj = w[v];
       Pythia8::Vec4 tem(wj.px(),wj.py(), wj.pz(), wj.e());*/
    
    Pythia8:: Vec4 final=operator+( temc,in);
    //tempdif=fabs(m_topmass - final.mCalc());
    tempdif=fabs(m_topmass - final.mCalc());
    
    if(tempdif < tpmass){
      tpmass=tempdif;
      bpos=a;
    } 
  }
  
  // matchedjets.push_back(w[wpos]);
  //returns best bjet followed by best w 
  // return btaggedjets[int(bpos)];
  return bpos;
}
