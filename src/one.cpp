#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TMath.h>

void one(){

   //TFile *file = new TFile("3tree.root","recreate");
   //TTree *tree = new TTree("tree","tree");

   Double_t x, y, z, dx, dy, dz;
   //tree->Branch("x", &x, "x/D");
   //tree->Branch("y", &y, "y/D");
   //tree->Branch("z", &z, "z/D");
   //tree->Branch("dx", &dx, "dx/D");
   //tree->Branch("dy", &dy, "dy/D");
   //tree->Branch("dz", &dz, "dz/D");
   Int_t N=0;

   Double_t x_th = 205;
   Double_t y_th = -110;
   TFile *f1 = new TFile("one.root", "recreate");
   TTree *tree1 = new TTree("tree1","tree1");
   tree1->Branch("x", &x, "x/D");
   tree1->Branch("y", &y, "y/D");
   tree1->Branch("z", &z, "z/D");
   tree1->Branch("dx", &dx, "dx/D");
   tree1->Branch("dy", &dy, "dy/D");
   tree1->Branch("dz", &dz, "dz/D");
   TGraph2D *g2d = new TGraph2D();
   g2d->SetNameTitle("g2d","Test Scan Counter; x [mm]; y [mm]; z [mm]");

   ifstream data("20150813.dat");
   while(data >> x >> y >> z >> dx >> dy >> dz){
      //tree->Fill();
      //N += 1;
      if(x > x_th) continue;
      if(y < y_th) continue;
      else{tree1->Fill();
      g2d->SetPoint(N++,x,y,z);
      }
   }
   data.close();

   //tree->Write();
   g2d->SetMarkerStyle(20);
   g2d->SetMarkerSize(0.02);
   g2d->SetMarkerColor(kBlack);
   g2d->Write();

   tree1->Write();

   //Int_t nevent = tree->GetEntries();
   Int_t nevent1 = tree1->GetEntries();
   //std::cout << "number of events : " << nevent << std::endl;
   std::cout << "number of events for one: " << nevent1 << std::endl;
   //std::cout << "N = " << N << std::endl;
}
