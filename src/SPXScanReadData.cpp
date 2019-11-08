#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TMath.h>

void SPXScanReadData() {

   TFile *file = new TFile("SPXScanReadData_190318.root","recreate");
   TTree *tree = new TTree("tree","tree");

   Double_t x, y, z;
   tree->Branch("x", &x, "x/D");
   tree->Branch("y", &y, "y/D");
   tree->Branch("z", &z, "z/D");
   
   Double_t a, b, c;

   ifstream data("190318_MEGII_TC_US_Final.asc");
   while(data >> x >> y >> z >> a >> b >> c) tree->Fill();
   data.close();

   tree->Write();

   Int_t nevent = tree->GetEntries();
   std::cout << "number of events : " << nevent <<std::endl;

}
