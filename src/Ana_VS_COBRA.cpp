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

void Ana_VS_COBRA()
{
   vector<TVector3> v_geo;
   vector<Double_t> phi_geo, theta_geo, psi_geo;
   const Int_t bufsize = 256;

   Int_t id,idx,order,hitorder;
   Double_t x,y,z,phi,theta,psi;
   TVector3 v;
   string str;
   ifstream ifs("data/Geometry.csv");
   if(!ifs){
      cout << "Error: Geometry.csv cannot open! Are you now at master directory?" << endl;
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
      cout << "Error: FARO_COBRA_calib2019.csv cannot open! Are you now at master directory?" << endl;
   } 
   while(getline(ifs,str)){
      sscanf(str.data(),"%lf,%lf",&x,&y);
      par.push_back(x);
      par_err.push_back(y);
   }
   ifs.close();
   const TVector3 v_par(par[0],par[1],par[2]);
 
   vector<TVector3> cv;
   ifs.open("outputfiles/Analysis.csv");
   if(!ifs){
      cout << "Error: File cannot open!" << endl;
   }
   while(getline(ifs,str)){
      sscanf(str.data(),"%lf,%lf,%lf",&x,&y,&z);
      v = TVector3(x,y,z);
      cv.push_back(v);
   }
   ifs.close();

   TCanvas *ch = new TCanvas("ch","ch",0,0,1000,800);
   TGraph2D *g2dr = new TGraph2D();
   g2dr->SetNameTitle("g2dr","z - #phi - dr plot; #phi[deg]; z[mm]; dr[mm]");
   g2dr->SetMarkerStyle(8);
   g2dr->SetMarkerSize(1.0);
   g2dr->SetMarkerColor(kGray+3);
   TGraph2D *g2dx = new TGraph2D();
   g2dx->SetNameTitle("g2dx","z - #phi - dx plot; #phi[deg]; z[mm]; dx[mm]");
   g2dx->SetMarkerStyle(8);
   g2dx->SetMarkerSize(1.0);
   g2dx->SetMarkerColor(kRed);
   TGraph2D *g2dy = new TGraph2D();
   g2dy->SetNameTitle("g2dy","z - #phi - dy plot; #phi[deg]; z[mm]; dy[mm]");
   g2dy->SetMarkerStyle(8);
   g2dy->SetMarkerSize(1.0);
   g2dy->SetMarkerColor(kYellow+1);
   TGraph2D *g2dz = new TGraph2D();
   g2dz->SetNameTitle("g2dz","z - #phi - dz plot; #phi[deg]; z[mm]; dz[mm]");
   g2dz->SetMarkerStyle(8);
   g2dz->SetMarkerSize(1.0);
   g2dz->SetMarkerColor(kBlue);
   TVector3 fv;
   for(Int_t i=0; i<48; i++){
      fv = cv[i];
      if(i <= 15){
         v = v_geo[511-16*i];
      }else if(i <= 31){
         v = v_geo[511-16*(i-16)-1];
      }else if(i <= 47){
         v = v_geo[511-16*(i-32)-15];
      }
      x = v.X();
      y = v.Y();
      z = v.Z();
      phi = TMath::ATan(y/x)*TMath::RadToDeg();
      fv = fv - v_par;
      fv.RotateX(par[3]*TMath::DegToRad());
      fv.RotateY(par[4]*TMath::DegToRad());
      fv.RotateZ(par[5]*TMath::DegToRad());
      fv = fv - v;
      g2dx->SetPoint(i,phi,z,fv.X());
      g2dy->SetPoint(i,phi,z,fv.Y());
      g2dz->SetPoint(i,phi,z,fv.Z());
      g2dr->SetPoint(i,phi,z,fv.Mag());
   }
   Double_t zmax = g2dx->GetYmax()+200;
   Double_t zmin = g2dx->GetYmin()-200;
   Double_t phimax = g2dx->GetXmax()+10;
   Double_t phimin = g2dx->GetXmin()-10;
   TPolyMarker3D *plane = new TPolyMarker3D();
   TRandom *rnd = new TRandom();
   rnd->SetSeed();
   for(Int_t i=0; i<10000; i++){
      y = zmin + (zmax-zmin) * rnd->Rndm();
      x = phimin + (phimax-phimin) * rnd->Rndm();
      z = 0.;
      plane->SetPoint(i,x,y,z);
   }
   plane->SetMarkerStyle(8);
   plane->SetMarkerSize(0.1);
   plane->SetMarkerColorAlpha(kGray+1,0.35);
   ch->Divide(2,2);
   ch->cd(1);
   g2dx->Draw("P");
   plane->Draw("same");
   g2dx->SetMargin(0.3);
   ch->cd(2);
   g2dy->Draw("P");
   plane->Draw("same");
   g2dy->SetMargin(0.3);
   g2dy->SetMaximum(10);
   g2dy->SetMinimum(-1.);
   ch->cd(3);
   g2dz->Draw("P");
   plane->Draw("same");
   g2dz->SetMargin(0.2);
   ch->cd(4);
   g2dr->Draw("P");
   plane->Draw("same");
   g2dr->SetMargin(0.2);
   g2dr->SetMaximum(10);
   g2dr->SetMinimum(-1.);
}

