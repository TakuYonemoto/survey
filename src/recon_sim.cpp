#include "TROOT.h"                                                                                  
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include <sstream>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TVector3.h>
#include <TMath.h>

void recon_sim()
{
   TFile *fr = new TFile("results_simu_kaizou.root");
   TTree *tree = (TTree*)fr->Get("tree");

   Double_t cx,cy,cz,ta[1000],tb[1000],tc[1000];

   tree->SetBranchAddress("cx",&cx);
   tree->SetBranchAddress("cy",&cy);
   tree->SetBranchAddress("cz",&cz);
   //tree->SetBranchAddress("ta",&ta);
   //tree->SetBranchAddress("tb",&tb);
   //tree->SetBranchAddress("tc",&tc);

   TGraph2D *g2drc[1000];
   TGraph2D *g2dor[1000];
   TGraph2D *gg = new TGraph2D();
   Double_t x,y,z;
   TVector3 v,yy,zz;
   Long64_t n;
   char* filename;
   filename = new char[50];
   TGraph2D *g2d;
   TFile *fs;
   Int_t nn = 3;

   for(Int_t k=0; k<nn; k++){
      fr->cd();
      tree->SetBranchAddress("ta",&ta[k]);
      tree->SetBranchAddress("tb",&tb[k]);
      tree->SetBranchAddress("tc",&tc[k]);
      tree->GetEntry(k);

      sprintf(filename, "./simu/sampleRR%04d.root",k);
      fs = new TFile(filename);
      g2d = (TGraph2D*)fs->Get("Graph2D");
      n = g2d->GetN();
      g2dor[k] = g2d;

      yy = TVector3(0,1,0);
      zz = TVector3(0,0,1);
      yy.RotateZ(ta[k] * TMath::DegToRad());
      zz.Rotate(tb[k] * TMath::DegToRad(), yy);

      for(Int_t i=0; i<n; i++){
         g2d->GetPoint(i,x,y,z);
         x = x - cx;
         y = y - cy;
         z = z - cz;
         v = TVector3(x,y,z);
         v.Rotate(-tc[k] * TMath::DegToRad(), zz);
         v.Rotate(-tb[k] * TMath::DegToRad(), yy);
         v.RotateZ(-ta[k] * TMath::DegToRad());
         x = v.X();
         y = v.Y();
         z = v.Z();
         gg->SetPoint(i,x,y,z);
      }
      g2drc[k] = gg;
   }
   fr->Close();

   TCanvas *c = new TCanvas("c","Graph2D",0,0,1000,800);
   c->Divide(2,3);
   for(Int_t k=0; k<nn*2; k++){
      c->cd(k+1);
      if(k % 2 == 0){
         g2dor[k/2]->Draw();
         g2dor[k/2]->SetMarkerStyle(20);
         g2dor[k/2]->SetMarkerSize(0.3);
      }
      else if(k % 2 == 1){
         g2drc[k/2]->Draw();
         g2drc[k/2]->SetMarkerStyle(20);
         g2drc[k/2]->SetMarkerSize(0.3);
      }
   }
}
