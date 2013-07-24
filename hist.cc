// Stdlib header file for input and output.
#include <iostream>

// Header file to access Pythia 8 program elements.
#include "Pythia.h"

// ROOT, for histogramming.
#include "TH1.h"

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT, for saving file.
#include "TFile.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "MyEvent.h"

//fastjet 
#include "fjClustering.h"

//Jets
#include "JetSeparation.h"
#include "MyTopEvent.h"
using namespace Pythia8; 

//functions 
bool isBhadron(int);
//void HFill(const Particle, std::string);
bool isElectron(const Particle &ptcl) { return ptcl.idAbs()==11;}
bool isTau(const Particle &ptcl) { return ptcl.idAbs() ==15;}
bool isMu(const Particle &ptcl) { return ptcl.idAbs()==13;}
bool isNu(const Particle &ptcl) { return ptcl.idAbs()==12 || ptcl.idAbs()==14 || ptcl.idAbs()==16; }
void fatal(TString msg) { printf("ERROR:\n\n  %s\n\n",msg.Data()); abort(); }
void PrintPtcl(const Particle &ptcl, TString comment="") {
  printf("  (pT,eta,phi,m) = (%6.1f GeV,%5.2f,%5.2f,%5.1f GeV) pdgID %4d : %s\n",
	 ptcl.pT(),ptcl.eta(),ptcl.phi(),ptcl.m(),ptcl.id(),comment.Data());
} 


//main function
int  main(int argc, char* argv[]) {
  
  bool partonMode=false;
  bool debug= true;

  // Create the ROOT application environment. 
  TApplication theApp("hist", &argc, argv);

  //generator
  Pythia pythia;
  
  //read settings from command file
  // if(partonMode){pythia.readFile("partonMode.cmnd");} else { pythia.readFile("particleMode.cmnd");}
  
  pythia.readString("SoftQCD:minBias = off");
  pythia.readString("SoftQCD:singleDiffractive = off");
  pythia.readString("SoftQCD:doubleDiffractive = off");
  pythia.readString("HardQCD:all = off");
  pythia.readString("PartonLevel:FSR = off");
  pythia.readString("PartonLevel:ISR = on ");
  pythia.readString("HadronLevel:all = on");
  pythia.readString("Top:qqbar2ttbar = on");
  pythia.readString("Top:gg2ttbar = on");
  
  if (partonMode) {
    pythia.readString("SoftQCD:minBias = on ");
    pythia.readString("SoftQCD:singleDiffractive = off");
    pythia.readString("SoftQCD:doubleDiffractive = off");
    pythia.readString("HardQCD:all = off");
    pythia.readString("PartonLevel:FSR = on");
    pythia.readString("PartonLevel:ISR = on ");
    pythia.readString("HadronLevel:all = on ");
    pythia.readString("Top:qqbar2ttbar = on");
    pythia.readString("Top:gg2ttbar = on");
  }
  
  //beam initialization (proton,proton, 7000 GeV)
  pythia.init(2212,2212,7000);
  
  //Create file on which histograms will be saved.
  TString of_name = partonMode ? "ttbar_partonLevel_histos.root" : "ttbar_hadronLevel_histograms.root";
  TFile* outFile = new TFile(of_name, "RECREATE");

 //Book histograms 
  TH1F* top_mass = new TH1F("m_top"," top mass;mass [GeV] ",100,140.,190.);
  TH1F* top_reconmass = new TH1F("reconmass_top","reconstructed mass of top", 100,0.,200.);
  TH1F* top_pt = new TH1F("pT_top","top quark pT;#it{p}_{T} [GeV];Frequency", 100, 0.  , 500.);
  TH1F* top_eta = new TH1F("eta_top","pseudrapitity of t;#eta^{top}", 100, -5., 5.);
  TH1F* top_y   = new TH1F("y_top",";#it{y}^{top}", 100, -5., 5.);

  TH1F* tbar_pt = new TH1F("pT_tbar","tbar pT", 100, -.5, 500.0);  
  TH1F* tbar_eta = new TH1F("eta_tbar","tbar pT", 100, -.5, 500.0);
  TH1F* tbar_y = new TH1F("y_tbar","bottom pT;it{p}_{T} [GeV]",100,0,100);
  TH1F* tbar_reconmass = new TH1F("reconmass_tbar","recon tbar pT;#it{p}_{T} [GeV]",100,0,100);
  TH1F* tbar_mass = new TH1F("mtbar"," tbar mass; mass [GeV]", 100,140.,190.);

  TH1F* whad = new TH1F("whadronic"," whad mass; mass [GeV]", 100,100.,190.);

 
  fjClustering* jetclustering   = new fjClustering(fastjet::antikt_algorithm, .25, fastjet::E_scheme, fastjet::Best); 
  JetMatching* jetcuts  = new JetMatching();
  MyTopEvent* mytop = new MyTopEvent();

  //number of events
  int nEv=1000; 
    
  /* struct mystruct{
    
    vector<Pythia8::Particle> particledata;
    vector<fastjet::PseudoJet> pseudojetdata;
    } data;*/

  // begin event loop 
  for (int iEv = 0; iEv < nEv; ++iEv) { 
   
    jetclustering->ClearJets();
    if (!pythia.next()) continue;//generate events.skip if necessary
    if (debug && iEv==0) {pythia.info.list(); pythia.event.list();} 
  
    // partons: the parton top and it's decay prodcuts (W and b)
    vector< Particle> tops, Ws, bs;
    
    // storage for B Hadrons
    vector<Particle> Bhadrons;
    
    // storage for neutrinos, leptons, and _all_ final particles (parton or hadron level particles depending on mode)
    vector<Particle>  all_stbl_ptcls, electrons, mus, taus, nus;    
    
    //counter for number of tops and positions in event record
    int ntop(0);
    vector<int> top_indices, windices;
    
    Event &event = pythia.event;
    
    
    /*
     *     Find final top quarks   
     *     
     */
    
  
    //begin top search
    for (size_t iptcl(0); iptcl < event.size(); ++iptcl) {    
      const Particle& particle = event[iptcl];

      int d1,d2;
      if(particle.idAbs()==6){
	
	// tops that  do not radiate a photon or "wiggle"
	d1=particle.daughter1(); d2=particle.daughter2();  
	if(d1==d2) continue;
	if(event[d1].idAbs()==6 || event[d2].idAbs()==6) continue;
	
	// Getting here means that we have found the right top.  
	top_indices.push_back(iptcl);
	ntop++;

	if (ntop > 2) {
	  pythia.event.list();
	  for (size_t j(0); j < top_indices.size();++j) printf("top %d has index %d\n",j,top_indices[j]);
	  fatal(Form("Event %d has %d good tops!?",iEv,ntop));
	}
      }      
    } //end top loop

   
    /*
     *    Find W decay products
     *        
     */

    
    for (int itop=0; itop < top_indices.size(); ++itop) {
     
      //position of tops
      int top_pos = top_indices[itop];
      const Particle& top = event[top_pos];
      
      // We have a top! it should go to a W and a b (most of the time)
      int W_index = top.daughter1(), b_index = top.daughter2();
      if (event[W_index].idAbs()!= 24) { W_index = top.daughter2(); b_index = top.daughter1(); } 
      const Particle &W = event[W_index];
      const Particle &b = event[b_index];

      if(b.idAbs()!=5) { printf(" We have top -> %d %d\n",W.id(),b.id()); continue;}
      if (W.idAbs()!=24 ) {printf("My W is not a W! its a: %d\n",W.id()); continue;}

      //store particles
      Ws.push_back(W);
      bs.push_back(b);
      
      //fill histograms
      Vec4 temp(W.p() + b.p());
      if(top.id() == 6){
	top_mass->Fill(temp.mCalc());
	top_eta->Fill(temp.eta());
	top_y->Fill(temp.rap());
      }
      else 
	tbar_mass->Fill(temp.mCalc());
      
      
      if(debug){ cout <<"daughters of top are: "<<event[W_index].id() <<" and  " <<event[b_index].id()<<"\n"<<endl;
	cout<<"w size() = "<<Ws.size()<<endl;}
      
    }
    

    /*
     *       Find final particles and store them in vectors
     *
     */

    
    for( size_t ithpt(0); ithpt < event.size(); ++ithpt){
      const Particle &particle = event[ithpt];
      const Particle &m1 = event[particle.mother1()];
      const Particle &m2= event[particle.mother2()];
      const Particle &m3= event[m1.mother1()];
      const Particle &m4 = event[m2.mother2()];

      if(! particle.isFinal()) continue;
     
      if ( particle.isLepton() && (m1.idAbs()==24 || m2.idAbs() ==24) &&  (m3.idAbs()==6 || m4.idAbs()==6)) {

	if( isNu(particle) & (m1.id()==24 || m2.id()==24))
	  nus.push_back(particle);
	
	if( isMu(particle) & (m1.id()==24 || m2.id()==24))
	   mus.push_back(particle);
	
	if( isElectron(particle) & (m1.id()==24 || m2.id()==24))
	  electrons.push_back(particle);

	continue;
      }
	
      if(isTau(particle)  && (m1.idAbs()==24 || m2.idAbs() ==24) && (m3.idAbs()==6 || m4.idAbs()==6)) {
	taus.push_back(particle);}

      if ((isNu(particle) || isMu(particle) || isElectron(particle)) && ( isTau(m1) || isTau(m2) ) ) continue; //leptons from tau decay

      if (isBhadron(particle.id())) 
	Bhadrons.push_back(particle);     
      
      jetclustering->push_back(particle);
    }

    if( electrons.size() ==0 && nus.size()==0 && taus.size()==0 && mus.size() == 0) continue; //skip this event
    if(nus.size()==0) continue;
    cout<<" eees"<<electrons.size()<<endl;
    cout<<" nuuu"<<nus.size()<<endl;
    cout<<" tuuu"<<taus.size()<<endl;
    cout<<" muuu"<<mus.size()<<endl;
    


    /*
     *      loop over jets and sort them
     *       
     */
     
     
     //cluster jets
     jetclustering->doClustering();
     
     if (debug  && iEv<5) {
       printf("\nEvent %d:\n",iEv);
       printf("  We found %d b-quarks\n",int(bs.size()));
       for (size_t i=0;i < bs.size();++i) PrintPtcl(bs[i],Form("b-quark %d",i));
       jetclustering->PrintJets();
     }
     
     
     /*
      *      Match jets to particles
      *      
      */

    
     vector<fastjet::PseudoJet> all_jets = jetclustering->GetJets();
     vector<fastjet::PseudoJet> btags,lightjets;
     
     jetcuts->SetParam(.5,2.4,3.0);
     all_jets= jetcuts->OverlapRemoval(mus,all_jets);
     
     jetcuts->SetParam(.5, 2.4, 3.0);
     all_jets= jetcuts->OverlapRemoval(electrons, all_jets);

     jetcuts->SetParam(.4, 2.4, 10.0); // delta R, etamax, ptmin
     btags = jetcuts->Match(bs, all_jets);
     lightjets = jetcuts->RemoveSubset(btags, all_jets);
     
     if(debug) {
       cout<<" \n size of all jets is"<< all_jets.size() <<endl;
       cout<<"\n number of matched bjets is: " <<btags.size()<< "\n"<<endl;
       for( int n(0) ; n < btags.size(); ++n) 
	 printf(" (pT,eta,phi,m) = (%6.1f GeV,%5.2f,%5.2f,%5.1f GeV) : \n",
		btags[n].pt(),
		btags[n].eta(),
		btags[n].phi(),
		btags[n].m());	     

       cout<<"\n lightjets are: \n"<<endl;
       for(size_t m(0); m < lightjets.size(); ++m)
	 printf(" (pT,eta,phi,m) = (%6.1f GeV,%5.2f,%5.2f,%5.1f GeV) : \n",
		lightjets[m].pt(),
		lightjets[m].eta(),
		lightjets[m].phi(),
		lightjets[m].m());
     }

  
     /*
      *     PseudoTop  
      *         
      */ 

    
     //send in btagged jets and lightjets(sorted)
     // Pythia8::Vec4 nuv = mytop->LeptonicW(nus, taus,mus, electrons);
     // double qun= mytop->Recon_Mass1(btags, lightjets, nuv, nus[0]);
     //if(qun<0) continue;
     //  whad->Fill( mytop->ReturnWmass());
     // top_reconmass->Fill(qun);

	    
  } //end of event loop
  //
  // Statistical summary
  if(debug) pythia.stat(); 
  
  //write to file and close it
  outFile->Write();
  outFile->Close();
  
  return 0;
  
} //end main function 

//check for b hadrons
bool isBhadron(int pdg) {
  if (pdg<0) pdg *=-1;
  if (pdg<500) return false;
  if (pdg/100%5 == 0 || pdg/1000==5) return true; 
  return false;
}

//void Hfill(const Particle particle, std::string name) {
 
  /*  if(std::string.compare("top")){
    top_mass->fill(particle.m());
    top_eta->fill(particle.eta());
    top_y->fill(particle.y()); 
  }
  else if(std::string.compare("tbar")){
    tbar_eta->fill(particle.eta());
    tbar_y->fill(particle.y());
    tbar_mass->fill(particle.m());
  }
  // else if ("toprecon")
  // top_reconmass->fill();
  // case "tbarrecon":
  // tbar_reconmass->fill();
  */
  //}
