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

void Initial_value()
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
 
   TCanvas *ch = new TCanvas("ch","ch",0,0,1000,800);
   TH1D *th1d[4];
   th1d[0] = new TH1D("de_x","difference_x",30,-10,10);
   th1d[1] = new TH1D("de_y","difference_y",30,-10,10);
   th1d[2] = new TH1D("de_z","difference_z",30,-10,10);
   th1d[3] = new TH1D("de_r","distance",15,0,10);
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
   for(Int_t i=0; i<16; i++){
      cv[i] = cv[i] - v_par;
      cv[i].RotateX(par[3]*TMath::DegToRad());
      cv[i].RotateY(par[4]*TMath::DegToRad());
      cv[i].RotateZ(par[5]*TMath::DegToRad());
      cv[i] = cv[i] - v_geo[511-16*i];
      cout << cv[i].X() << ", " << cv[i].Y() << ", " << cv[i].Z() << endl;
      th1d[0]->Fill(cv[i].X());
      th1d[1]->Fill(cv[i].Y());
      th1d[2]->Fill(cv[i].Z());
      th1d[3]->Fill(cv[i].Mag());
   }
   for(Int_t i=16; i<32; i++){
      cv[i] = cv[i] - v_par;
      cv[i].RotateX(par[3]*TMath::DegToRad());
      cv[i].RotateY(par[4]*TMath::DegToRad());
      cv[i].RotateZ(par[5]*TMath::DegToRad());
      cv[i] = cv[i] - v_geo[511-16*(i-16)-1];
      cout << cv[i].X() << ", " << cv[i].Y() << ", " << cv[i].Z() << endl;
      th1d[0]->Fill(cv[i].X());
      th1d[1]->Fill(cv[i].Y());
      th1d[2]->Fill(cv[i].Z());
      th1d[3]->Fill(cv[i].Mag());
   }
   for(Int_t i=32; i<48; i++){
      cv[i] = cv[i] - v_par;
      cv[i].RotateX(par[3]*TMath::DegToRad());
      cv[i].RotateY(par[4]*TMath::DegToRad());
      cv[i].RotateZ(par[5]*TMath::DegToRad());
      cv[i] = cv[i] - v_geo[511-16*(i-32)-15];
      cout << cv[i].X() << ", " << cv[i].Y() << ", " << cv[i].Z() << endl;
      th1d[0]->Fill(cv[i].X());
      th1d[1]->Fill(cv[i].Y());
      th1d[2]->Fill(cv[i].Z());
      th1d[3]->Fill(cv[i].Mag());
   }
   ch->Divide(2,2);
   ch->cd(1);
   th1d[0]->Draw();
   ch->cd(2);
   th1d[1]->Draw();
   ch->cd(3);
   th1d[2]->Draw();
   ch->cd(4);
   th1d[3]->Draw();



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
      y = -3.1 + (generator->Rndm()*6.2);
      z = -60. + (generator->Rndm()*120.);
      pm3->SetPoint(N++,x,y,z);
   }
   for(Int_t i=0; i < 3000; i++){
      x = -20. + (generator->Rndm()*40.);
      y = -3.1 + (generator->Rndm()*6.2);
      z = 0;
      pm3->SetPoint(N++,x,y,z);
   }
   for(Int_t i=0; i < 3000; i++){
      x = -20. + (generator->Rndm()*40.);
      y = 0;
      z = -60. + (generator->Rndm()*120.);
      pm3->SetPoint(N++,x,y,z);
   }
   TPolyMarker3D *pm3h = new TPolyMarker3D();
   pm3h->SetName("pm3_origin_5cm");
   pm3h->SetMarkerStyle(8);
   pm3h->SetMarkerSize(0.1);
   pm3h->SetMarkerColor(kBlue);
   N = 0;
   for(Int_t i=0; i < 3000; i++){
      x = 0;
      y = -3.1 + (generator->Rndm()*6.2);
      z = -60. + (generator->Rndm()*120.);
      pm3h->SetPoint(N++,x,y,z);
   }
   for(Int_t i=0; i < 4000; i++){
      x = -25. + (generator->Rndm()*50.);
      y = -3.1 + (generator->Rndm()*6.2);
      z = 0;
      pm3h->SetPoint(N++,x,y,z);
   }
   for(Int_t i=0; i < 4000; i++){
      x = -25. + (generator->Rndm()*50.);
      y = 0;
      z = -60. + (generator->Rndm()*120.);
      pm3h->SetPoint(N++,x,y,z);
   }

   TFile *inputfile = new TFile("/Users/taku/Downloads/survey_data/kiridashi.root","read");
   TFile *outputfile = new TFile("outputfiles/Initial_value.root","recreate");
   pm3->Write();
   pm3h->Write();
   Double_t y_mean = 0.;
   for(Int_t iPixel = 0; iPixel<256; iPixel++){
      if(iPixel%20) continue;
      y_mean = 0.;
      cout << iPixel << endl;
      TTree *tree_copy = ((TTree*)inputfile->Get(Form("t3_%d",iPixel)))->CloneTree();
      if(!tree_copy){cout << "input tree cannot be read." << endl;}
      Int_t nEntry = tree_copy->GetEntries();
      tree_copy->SetBranchAddress("x",&x);
      tree_copy->SetBranchAddress("y",&y);
      tree_copy->SetBranchAddress("z",&z);
      tree_copy->SetBranchStatus("z2",0);
      tree_copy->SetBranchStatus("phi",0);
      tree_copy->SetBranchStatus("phi2",0);
      Int_t USPixel = 511 - iPixel;
      TTree *transtree = new TTree(Form("trans_%d",USPixel),Form("trans_%d",USPixel));
      transtree->Branch("x",&x,"x/D");
      transtree->Branch("y",&y,"y/D");
      transtree->Branch("z",&z,"z/D");
      TPolyMarker3D *pm3x = new TPolyMarker3D();
      pm3x->SetName(Form("pm3x_%d",USPixel));
      pm3x->SetMarkerStyle(8);
      pm3x->SetMarkerSize(0.1);
      pm3x->SetMarkerColor(kRed);
      N = 0;
      for(Int_t iEntry=0; iEntry < nEntry; iEntry++){
         tree_copy->GetEntry(iEntry);
         TVector3 v_FARO(x,y,z);
         TVector3 v_Rot(v_FARO - v_par);
         v_Rot.RotateX(par[3]*TMath::DegToRad());
         v_Rot.RotateY(par[4]*TMath::DegToRad());
         v_Rot.RotateZ(par[5]*TMath::DegToRad());
         TVector3 v_Rot2(v_Rot - v_geo[USPixel]);
         v_Rot2.RotateZ(-psi_geo[USPixel]);
         v_Rot2.RotateX(-theta_geo[USPixel]);
         v_Rot2.RotateZ(-phi_geo[USPixel]);
         x = v_Rot2.X();
         y = v_Rot2.Y();
         z = v_Rot2.Z();
         transtree->Fill();
         if(x < -19 && -22 < x){
            if((USPixel >= 288 && USPixel <= 302) && USPixel%2 == 0) continue;
            else if((USPixel >= 384 && USPixel <= 399) && USPixel%2 != 0) continue;
            else if(USPixel >= 304 && USPixel <= 383) continue;
            pm3x->SetPoint(N++,x,y,z);
            y_mean += y;
         }
         if(x < -24 && -27 < x){
            if((USPixel >= 288 && USPixel <= 302) && USPixel%2 != 0) continue;
            else if((USPixel >= 384 && USPixel <= 399) && USPixel%2 == 0) continue;
            else if(USPixel <= 304 && USPixel >= 383) continue;
            pm3x->SetPoint(N++,x,y,z);
            y_mean += y;
         }
      }
      transtree->Write();
      if(N!= 0){
         y_mean /= N;
         TPolyMarker3D *pm3y = new TPolyMarker3D();
         pm3y->SetName(Form("pm3y_%d",USPixel));
         pm3y->SetMarkerStyle(8);
         pm3y->SetMarkerSize(0.2);
         pm3y->SetMarkerColor(kGreen+2);
         for(Int_t iPoint=0; iPoint<10000; iPoint++){
            x = -20. + (generator->Rndm()*40.);
            y = y_mean;
            z = -60. + (generator->Rndm()*120.);
            if((USPixel >= 288 && USPixel <= 302) && USPixel%2 == 0) x -= 5.;
            else if((USPixel >= 384 && USPixel <= 399) && USPixel%2 != 0) x -= 5.;
            else if(USPixel >= 304 && USPixel <= 383) x -= 5;
            pm3y->SetPoint(iPoint,x,y,z);
         }
         pm3x->Write();
         pm3y->Write();
      }
   }
   outputfile->Close();
}

