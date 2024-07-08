#include <TROOT.h>
#include <TMath.h>
#include <TVector3.h>
#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TGraph2D.h>
#include <TPolyMarker3D.h>
#include <TRandom.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;


void Autokiridashi_FARO_CAD()
{
   vector<Double_t> x_geo, y_geo, z_geo;
   vector<Double_t> phi_geo, theta_geo, psi_geo;
   vector<Double_t> par, par_err;
   Double_t x,y,z,phi,theta,psi;
   
   Int_t id,idx,order,hitorder;
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
      x_geo.push_back(x);
      y_geo.push_back(y);
      z_geo.push_back(z);
      phi_geo.push_back(phi);
      theta_geo.push_back(theta);
      psi_geo.push_back(psi);
   }
   ifs.close();

   ifs.open("outputfiles/FARO_CAD_calib.csv");
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
   pm3->SetName("pm3_origin");
   pm3->SetMarkerStyle(8);
   pm3->SetMarkerSize(0.1);
   pm3->SetMarkerColor(kBlue);
   Int_t N = 0;
   for(Int_t i=0; i < 3000; i++){
      x = 0;
      y = 0;
      z = -60. + (generator->Rndm()*120.);
      pm3->SetPoint(N++,x,y,z);
   }
   for(Int_t i=0; i < 2000; i++){
      x = 0;
      y = -6. + (generator->Rndm()*12.);
      z = 0;
      pm3->SetPoint(N++,x,y,z);
   }
   for(Int_t i=0; i < 3000; i++){
      x = -20. + (generator->Rndm()*40.);
      y = 0;
      z = 0;
      pm3->SetPoint(N++,x,y,z);
   }
   
   TVector3 v;
   Double_t xt,yt,zt;
   TFile *data = new TFile("/Users/taku/Downloads/survey_data/SPXScanReadData_190318_US.root","read");
   TFile *fout = new TFile("/Users/taku/Downloads/survey_data/Autokiridashi_FARO_CAD.root","recreate");
   TTree *td = ((TTree*)data->Get("tree"))->CloneTree();;
   td->SetBranchStatus("*",0); 
   td->SetBranchStatus("x",1); 
   td->SetBranchStatus("y",1);
   td->SetBranchStatus("z",1);
   td->SetBranchAddress("x",&x); 
   td->SetBranchAddress("y",&y);
   td->SetBranchAddress("z",&z);
   Int_t nEntry = td->GetEntries();
   pm3->Write();
   Double_t x_min,y_max,y_min;

   for(Int_t iPixel = 0; iPixel < 30; iPixel++){
      TGraph2D *g2d = new TGraph2D();
      N = 0;
      cout << iPixel << endl;
      Int_t USPixel = iPixel + 256;
      TTree *t3 = new TTree(Form("t3_%d",iPixel),Form("t3_%d",iPixel));
      t3->Branch("x",&x,"x/D");
      t3->Branch("y",&y,"y/D");
      t3->Branch("z",&z,"z/D");
      t3->Branch("xt",&xt,"xt/D");
      t3->Branch("yt",&yt,"yt/D");
      t3->Branch("zt",&zt,"zt/D");
      for(Int_t iEntry = 0; iEntry < nEntry; iEntry++){
         td->GetEntry(iEntry);
         TVector3 v_Rot(x,y,z);
         TVector3 v_geo(x_geo[USPixel],y_geo[USPixel],z_geo[USPixel]);
         v_Rot = v_Rot - v_par;
         v_Rot.RotateX(par[3]*TMath::DegToRad());
         v_Rot.RotateY(par[4]*TMath::DegToRad());
         v_Rot.RotateZ(par[5]*TMath::DegToRad());
         v_Rot = v_Rot - v_geo;
         v_Rot.RotateZ(-psi_geo[USPixel]);
         v_Rot.RotateX(-theta_geo[USPixel]);
         v_Rot.RotateZ(-phi_geo[USPixel]);
         xt = v_Rot.X();
         yt = v_Rot.Y();
         zt = v_Rot.Z();
         if(fabs(xt) < 27. && fabs(yt) < 10. && fabs(zt) < 80.){
            t3->Fill();
            if(xt < -19 && -21 < xt){
               g2d->SetPoint(N++,xt,yt,zt);
            }
         }
      }
      t3->Write();
      if(N!= 0){
         y_max = g2d->GetYmax();
         y_min = g2d->GetYmin();
         TPolyMarker3D *pm3y = new TPolyMarker3D();
         pm3y->SetName(Form("pm3_%d",iPixel));
         pm3y->SetMarkerStyle(8);
         pm3y->SetMarkerSize(0.1);
         pm3y->SetMarkerColor(kRed);
         for(Int_t iPoint=0; iPoint < N; iPoint++){
            g2d->GetPoint(iPoint,xt,yt,zt);
            pm3y->SetPoint(iPoint,xt,yt,zt);
         }
         for(Int_t iPoint=0; iPoint<10000; iPoint++){
            x = -20. + (generator->Rndm()*40.);
            y = (y_max + y_min)/2;
            z = -60. + (generator->Rndm()*120.);
            pm3y->SetPoint(iPoint+N,x,y,z);
         }
         pm3y->Write();
      }
   }
   fout->Close();
}

