#include "MyHisto.h"
#include "pythia.h"
#include "fjClustering.h"
#include "MyEvent.h"

//this needs to be more general
void Create( bool partonMode){

 //Create file on which histograms will be saved.
  TString of_name = partonMode ? "ttbar_partonLevel_histos.root" : "ttbar_hadronLevel_histograms.root";
  TFile* outFile = new TFile(of_name, "RECREATE");

 //Book histograms 
  TH1F* top_mass = new TH1F("m_top"," top mass",100,160.,190.);
  TH1F* top_reconst = new TH1F("reconmass_top","reconstructed mass of top", 100,150.,200.);
  TH1F* top_pt = new TH1F("pT_top","top quark pT;#it{p}_{T} [GeV];Frequency", 100, 0.  , 500.);
  TH1F* top_eta = new TH1F("eta_top","pseudrapitity of t;#eta^{top}", 100, -5., 5.);
  TH1F* top_y   = new TH1F("y_top",";#it{y}^{top}", 100, -5., 5.);

  TH1F* tbar_pt = new TH1F("pT_tbar","tbar pT", 100, -.5, 500.0);
  TH1F* bquarkPt = new TH1F("pT_bottom","bottom pT;it{p}_{T} [GeV]",100,0,100);
  TH1F* bbarquarkPt = new TH1F("pT_bbar","bbar pT;#it{p}_{T} [GeV]",100,0,100);

}

void HFill(const Pythia8::Particle &ptcl) {
  
  

}

void CloseFile(){
  //write to and close root file
  outFile->Write();
  outFile->Close();
}

