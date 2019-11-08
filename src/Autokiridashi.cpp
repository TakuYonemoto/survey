#include <TMath.h>
#include <TVector3.h>
#include <TTree.h>
#include <TFile.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

void Autokiridashi()
{   
   vector<Double_t> x_geo, y_geo, z_geo;
   vector<Double_t> phi_geo, theta_geo, psi_geo;

   Int_t id,idx,order,hitorder;
   Double_t x,y,z,phi,theta,psi;
   TVector3 v;
   string str;
   ifstream ifs("data/Geometry.csv");
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
   Double_t xt,yt,zt;

   TFile *data = new TFile("/Users/taku/Downloads/survey_data/SPXScanReadData_190318_US.root","read");
   TFile *fout = new TFile("/Users/taku/Downloads/survey_data/Autokiridashi.root","recreate");
   TTree *td = ((TTree*)data->Get("tree"))->CloneTree();;
   td->SetBranchStatus("*",0); 
   td->SetBranchStatus("x",1); 
   td->SetBranchStatus("y",1);
   td->SetBranchStatus("z",1);
   td->SetBranchAddress("x",&x); 
   td->SetBranchAddress("y",&y);
   td->SetBranchAddress("z",&z);
   Int_t nEntry = td->GetEntries();

   for(Int_t iPixel = 0; iPixel < 256; iPixel++){
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
         TVector3 v_FARO(x,y,z);
         TVector3 v_geo(x_geo[USPixel],y_geo[USPixel],z_geo[USPixel]);
         TVector3 v_Rot(v_FARO - v_par);
         v_Rot.RotateX(par[3]*TMath::DegToRad());
         v_Rot.RotateY(par[4]*TMath::DegToRad());
         v_Rot.RotateZ(par[5]*TMath::DegToRad());
         TVector3 v_Rot2(v_Rot - v_geo);
         v_Rot2.RotateZ(-psi_geo[USPixel]);
         v_Rot2.RotateX(-theta_geo[USPixel]);
         v_Rot2.RotateZ(-phi_geo[USPixel]);
         xt = v_Rot2.X();
         yt = v_Rot2.Y();
         zt = v_Rot2.Z();
         if(fabs(xt) < 27. && fabs(yt) < 7. && fabs(zt) < 70.){
            t3->Fill();
         }
      }
      t3->Write();
   }
   fout->Close();
}

