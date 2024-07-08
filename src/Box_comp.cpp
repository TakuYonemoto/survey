#include <TMath.h>
#include <TVector3.h>
#include <TTree.h>
#include <TFile.h>
#include <TPolyMarker3D.h>
#include <TGraph2D.h>
#include <TRandom.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

void Box_comp()
{
   vector<TVector3> v_geo;
   vector<Double_t> phi_geo, theta_geo, psi_geo;
   const Int_t bufsize = 256;

   Int_t id,idx,order,hitorder;
   Double_t x,y,z,phi,theta,psi;
   TVector3 v;
   string str;
   ifstream ifs("data/PixelPositions.csv");
   if(!ifs){
      cout << "Error: File cannot open!" << endl;
   } 
   while(getline(ifs,str)){
      sscanf(str.data(), "%d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%d,%d",&id,&idx,&x,&y,&z,&phi,&theta,&psi,&order,&hitorder);
      x *= 10.;
      y *= 10.;
      z *= 10.;
      v = TVector3(x,y,z);
      v_geo.push_back(v);
      phi_geo.push_back(phi);
      theta_geo.push_back(theta);
      psi_geo.push_back(psi);
   }
   ifs.close();

   vector<Double_t> par, par_err;
   ifs.open("outputfiles/FARO_COBRA_calib2019.csv");
   if(!ifs){
      cout << "Error: File cannot open!" << endl;
   } 
   while(getline(ifs,str)){
      sscanf(str.data(),"%lf,%lf",&x,&y);
      par.push_back(x);
      par_err.push_back(y);
   }
   ifs.close();
   const TVector3 v_par(par[0],par[1],par[2]);

   TRandom *generator = new TRandom();
   generator->SetSeed();              

   TPolyMarker3D *pm3 = new TPolyMarker3D();
   pm3->SetMarkerStyle(8);
   pm3->SetMarkerSize(0.1);
   pm3->SetMarkerColor(kRed);
   TFile *boxmake1 = new TFile("outputfiles/boxmake1.root", "read");
   TPolyMarker3D *ppm3 = (TPolyMarker3D*)boxmake1->Get("pm3");
   Int_t N = ppm3->GetN();
   for(Int_t i=0; i<N; i++){
      ppm3->GetPoint(i,x,y,z);
      pm3->SetPoint(i,y+8.154485,z,x);
   }
   boxmake1->Close();

   TFile *inputfile = new TFile("/Users/taku/Downloads/survey_data/kiridashi.root","read");
   TFile *outputfile = new TFile("/Users/taku/Downloads/survey_data/Box_comp.root","recreate");
   TPolyMarker3D *PM3[256];
   TTree *transtree[256], *tree_copy[256];

   TCanvas *c = new TCanvas("c","c",0,0,600,400);
   Int_t nEntries;
   TVector3 v_Rot;

   for(Int_t iPixel = 0; iPixel < 2; iPixel++){
      //if(iPixel%20) continue;
      cout << iPixel << endl;

      tree_copy[iPixel] = ((TTree*)inputfile->Get(Form("t3_%d",iPixel)))->CloneTree();
      if(!tree_copy[iPixel]){cout << "input tree cannot be read." << endl;}
      tree_copy[iPixel]->SetBranchAddress("x",&x);
      tree_copy[iPixel]->SetBranchAddress("y",&y);
      tree_copy[iPixel]->SetBranchAddress("z",&z);
      nEntries = tree_copy[iPixel]->GetEntries();

      transtree[iPixel] = new TTree(Form("trans_%d",iPixel),Form("trans_%d",iPixel));
      transtree[iPixel]->Branch("x",&x,"x/D");
      transtree[iPixel]->Branch("y",&y,"y/D");
      transtree[iPixel]->Branch("z",&z,"z/D");
      for(Int_t i=0; i < nEntries; i++){
         tree_copy[iPixel]->GetEntry(i);
         v_Rot = TVector3(x,y,z);
         v_Rot = v_Rot - v_par;
         v_Rot.RotateX(par[3]*TMath::DegToRad());
         v_Rot.RotateY(par[4]*TMath::DegToRad());
         v_Rot.RotateZ(par[5]*TMath::DegToRad());
         x = v_Rot.X();
         y = v_Rot.Y();
         z = v_Rot.Z();
         transtree[iPixel]->Fill();
      }
      transtree[iPixel]->Write();

      Int_t USPixel = 511 - iPixel;
      PM3[iPixel] = new TPolyMarker3D();
      PM3[iPixel]->SetName(Form("PM3_%d",USPixel));
      PM3[iPixel]->SetMarkerStyle(8);
      PM3[iPixel]->SetMarkerSize(0.3);
      PM3[iPixel]->SetMarkerColor(kRed);
      
      for(Int_t i=0; i < N; i++){
         pm3->GetPoint(i,x,y,z);
         v_Rot = TVector3(x,y,z);
         v_Rot.RotateZ(phi_geo[USPixel]);
         v_Rot.RotateX(theta_geo[USPixel]);
         v_Rot.RotateZ(psi_geo[USPixel]);
         v_Rot = v_Rot + v_geo[USPixel];
         x = v_Rot.X();
         y = v_Rot.Y();
         z = v_Rot.Z();
         PM3[iPixel]->SetPoint(i,x,y,z);
      }
      PM3[iPixel]->Write();
   }
   transtree[0]->Draw("z:y:x");
   PM3[0]->Draw("same");
   outputfile->Close();
}
