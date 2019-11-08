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

using namespace TMath;

/// fit-box size ///
const Double_t hx = 61.9;   //mm
const Double_t Hy0 = 56.35; //mm
const Double_t hz2 = 3.1;    //mm
const Double_t hz = 5.25;  //mm
const Double_t tanyz = (TMath::Tan(TMath::ASin((hz-hz2)/Hy0)));
const Double_t hy = (hz-hz2)/tanyz/2; //mm
const Double_t S1 = (2*hz2 + 2*hz) * hy / 2; //mm^2
const Double_t S2 = 2*hz2 * 2*hx;            //mm^2
const Double_t S3 = 2*hx * Hy0;            //mm^2
const Double_t S4 = 2*hx * 2*hz;             //mm^2
const Double_t S  = 2*S1 + 2*S3 + S2 + S4; 

const Double_t hX = 60.;  //mm
const Double_t hY = 20.;  //mm
const Double_t hZ = 2.5;  //mm

void prerotate()
{ 
   TFile *f = new TFile("/Users/taku/Downloads/survey_data/kiridashi.root");
   TTree *t0[16];
   TTree *t1[16];
   TTree *t15[16];
   Int_t t3num;
   for(Int_t i=0; i<16; i++){
      t3num = i * 16;
      t0[i] = (TTree*)f->Get(Form("t3_%d",t3num));
   }
   for(Int_t i=0; i<16; i++){
      t3num = i * 16 + 1;
      t1[i] = (TTree*)f->Get(Form("t3_%d",t3num));
   }
   for(Int_t i=0; i<16; i++){
      t3num = i * 16 + 15;
      t15[i] = (TTree*)f->Get(Form("t3_%d",t3num));
   }
   Long64_t n = 0;

   Double_t x,y,z;
   Double_t cx0[16],cy0[16],cz0[16],tx0[16],ty0[16],tz0[16];
   Double_t cx1[16],cy1[16],cz1[16],tx1[16],ty1[16],tz1[16];
   Double_t cx15[16],cy15[16],cz15[16],tx15[16],ty15[16],tz15[16];
   TFile *fo = new TFile("outputfiles/prerotate.root","recreate");
   TTree *tree0[16];
   for(Int_t i=0; i<16; i++){
      tree0[i] = t0[i]->CloneTree();
      tree0[i]->Write();
   }
   TTree *tree1[16];
   for(Int_t i=0; i<16; i++){
      tree1[i] = t1[i]->CloneTree();
      tree1[i]->Write();
   }
   TTree *tree15[16];
   for(Int_t i=0; i<16; i++){
      tree15[i] = t15[i]->CloneTree();
      tree15[i]->Write();
   }
   TTree *tree2 = new TTree("tree2","tree2");
   tree2->Branch("cx0", cx0, "cx0[16]/D"); 
   tree2->Branch("cy0", cy0, "cy0[16]/D"); 
   tree2->Branch("cz0", cz0, "cz0[16]/D");
   tree2->Branch("tx0", tx0, "tx0[16]/D"); 
   tree2->Branch("ty0", ty0, "ty0[16]/D"); 
   tree2->Branch("tz0", tz0, "tz0[16]/D");
   TTree *tree3 = new TTree("tree3","tree3");
   tree3->Branch("cx1", cx1, "cx1[16]/D"); 
   tree3->Branch("cy1", cy1, "cy1[16]/D"); 
   tree3->Branch("cz1", cz1, "cz1[16]/D");
   tree3->Branch("tx1", tx1, "tx1[16]/D"); 
   tree3->Branch("ty1", ty1, "ty1[16]/D"); 
   tree3->Branch("tz1", tz1, "tz1[16]/D");
   tree3->Branch("cx15", cx15, "cx15[16]/D"); 
   tree3->Branch("cy15", cy15, "cy15[16]/D"); 
   tree3->Branch("cz15", cz15, "cz15[16]/D");
   tree3->Branch("tx15", tx15, "tx15[16]/D"); 
   tree3->Branch("ty15", ty15, "ty15[16]/D"); 
   tree3->Branch("tz15", tz15, "tz15[16]/D");


   //t3_0
   cx0[0] = 52.7; cy0[0] = 51.5; cz0[0] = -53.5;
   tx0[0] = -45.; ty0[0] = 8.; tz0[0] = -82.;
   //t3_1
   cx1[0] = 69.3; cy1[0] = 111.2; cz1[0] = -23.5;
   tx1[0] = -43.5; ty1[0] = 13.; tz1[0] = -74.;
   //t3_2
   /*
   cx2[0] = 0.;
   cy2[0] = 0.;
   cz2[0] = 0.;
   tx2[0] = 0.;
   ty2[0] = 0.;
   tz2[0] = 0.;
   */
   //t3_15
   cx15[0] = 675.; cy15[0] = 65.; cz15[0] = -18.;
   tx15[0] = -136.; ty15[0] = -10.; tz15[0] = 100.;

   //t3_16
   cx0[1] = 54.; cy0[1] = 53.5; cz0[1] = -108.5;
   tx0[1] = -45.; ty0[1] = 8.; tz0[1] = -82.;
   //t3_17
   cx1[1] = 69.3; cy1[1] = 108.; cz1[1] = -79.;
   tx1[1] = -43.5; ty1[1] = 13.; tz1[1] = -74.5;
   //t3_31
   cx15[1] = 675.; cy15[1] = 65.; cz15[1] = -72.;
   tx15[1] = -135.; ty15[1] = -8.; tz15[1] = 100.;

   //t3_32
   cx0[2] = 55.; cy0[2] = 55.; cz0[2] = -162.;
   tx0[2] = -45.; ty0[2] = 8.; tz0[2] = -82.;
   //t3_33
   cx1[2] = 71.; cy1[2] = 110.; cz1[2] = -131.5;
   tx1[2] = -43.5; ty1[2] = 13.; tz1[2] = -74.5;
   //t3_47
   cx15[2] = 675.; cy15[2] = 65.; cz15[2] = -126.;
   tx15[2] = -135.; ty15[2] = -8.; tz15[2] = 100.;

   //t3_48
   cx0[3] = 55.; cy0[3] = 55.; cz0[3] = -216.;
   tx0[3] = -45.; ty0[3] = 8.; tz0[3] = -82.;
   //t3_49
   cx1[3] = 71.; cy1[3] = 110.; cz1[3] = -185.5;
   tx1[3] = -43.5; ty1[3] = 13.; tz1[3] = -74.5;
   //t3_63
   cx15[3] = 675.; cy15[3] = 65.; cz15[3] = -180.;
   tx15[3] = -135.; ty15[3] = -8.; tz15[3] = 100.;
   
   //t3_64
   cx0[4] = 55.; cy0[4] = 55.; cz0[4] = -272.;
   tx0[4] = -45.; ty0[4] = 8.; tz0[4] = -82.;
   //t3_65
   cx1[4] = 73.; cy1[4] = 110.; cz1[4] = -241.5;
   tx1[4] = -43.5; ty1[4] = 13.; tz1[4] = -74.;
   //t3_79
   cx15[4] = 675.; cy15[4] = 65.; cz15[4] = -237.;
   tx15[4] = -135.; ty15[4] = -8.; tz15[4] = 100.;

   //t3_80
   cx0[5] = 55.; cy0[5] = 55.; cz0[5] = -326.;
   tx0[5] = -45.; ty0[5] = 8.; tz0[5] = -82.;
   //t3_81
   cx1[5] = 73.; cy1[5] = 110.; cz1[5] = -299.;
   tx1[5] = -43.5; ty1[5] = 13.; tz1[5] = -74.;
   //t3_97
   cx15[5] = 676.; cy15[5] = 64.; cz15[5] = -292.;
   tx15[5] = -135.; ty15[5] = -8.; tz15[5] = 99.;

   //t3_96
   cx0[6] = 55.; cy0[6] = 55.; cz0[6] = -381.5;
   tx0[6] = -45.; ty0[6] = 8.; tz0[6] = -82.;
   //t3_97
   cx1[6] = 73.; cy1[6] = 110.; cz1[6] = -353.5;
   tx1[6] = -43.5; ty1[6] = 13.; tz1[6] = -74.;
   //t3_111
   cx15[6] = 676.; cy15[6] = 64.; cz15[6] = -347.;
   tx15[6] = -135.; ty15[6] = -8.; tz15[6] = 99.;

   //t3_112
   cx0[7] = 55.; cy0[7] = 55.; cz0[7] = -435.;
   tx0[7] = -45.; ty0[7] = 8.; tz0[7] = -82.;
   //t3_113
   cx1[7] = 73.; cy1[7] = 110.; cz1[7] = -408.;
   tx1[7] = -43.5; ty1[7] = 13.; tz1[7] = -74.;
   //t3_127
   cx15[7] = 676.; cy15[7] = 63.; cz15[7] = -401.;
   tx15[7] = -135.; ty15[7] = -8.; tz15[7] = 99.;

   //t3_128
   cx0[8] = 56.5; cy0[8] = 55.; cz0[8] = -492.;
   tx0[8] = -45.; ty0[8] = 8.; tz0[8] = -82.;
   //t3_129 
   cx1[8] = 73.; cy1[8] = 110.; cz1[8] = -462.;
   tx1[8] = -43.5; ty1[8] = 13.; tz1[8] = -74.;
   //t3_143
   cx15[8] = 678.; cy15[8] = 64.; cz15[8] = -456.;
   tx15[8] = -135.; ty15[8] = -8.; tz15[8] = 99.;

   //t3_144
   cx0[9] = 56.5; cy0[9] = 55.; cz0[9] = -546.;
   tx0[9] = -45.; ty0[9] = 8.; tz0[9] = -82.;
   //t3_145
   cx1[9] = 75.; cy1[9] = 110.; cz1[9] = -517.2;
   tx1[9] = -43.5; ty1[9] = 13.; tz1[9] = -74.;
   //t3_159
   cx15[9] = 678.; cy15[9] = 63.; cz15[9] = -510.;
   tx15[9] = -135.; ty15[9] = -8.; tz15[9] = 99.;

   //t3_160
   cx0[10] = 58.5; cy0[10] = 55.; cz0[10] = -602.;
   tx0[10] = -45.; ty0[10] = 8.; tz0[10] = -82.;
   //t3_161
   cx1[10] = 77.; cy1[10] = 110.; cz1[10] = -573.2;
   tx1[10] = -43.5; ty1[10] = 13.; tz1[10] = -74.;
   //t3_175
   cx15[10] = 680.; cy15[10] = 62.; cz15[10] = -566.;
   tx15[10] = -135.; ty15[10] = -8.; tz15[10] = 99.;

   //t3_176
   cx0[11] = 58.5; cy0[11] = 55.; cz0[11] = -656.;
   tx0[11] = -44.; ty0[11] = 8.; tz0[11] = -82.;
   //t3_177
   cx1[11] = 77.; cy1[11] = 110.; cz1[11] = -627.2;
   tx1[11] = -43.5; ty1[11] = 13.; tz1[11] = -74.;
   //t3_191
   cx15[11] = 680.; cy15[11] = 62.; cz15[11] = -618.;
   tx15[11] = -135.; ty15[11] = -8.; tz15[11] = 99.;

   //t3_192
   cx0[12] = 61.5; cy0[12] = 55.; cz0[12] = -713.;
   tx0[12] = -44.; ty0[12] = 8.; tz0[12] = -82.;
   //t3_193
   cx1[12] = 77.; cy1[12] = 110.; cz1[12] = -682.2;
   tx1[12] = -43.5; ty1[12] = 13.; tz1[12] = -74.;
   //t3_207
   cx15[12] = 681.; cy15[12] = 62.; cz15[12] = -674.;
   tx15[12] = -135.; ty15[12] = -8.; tz15[12] = 99.;

   //t3_208
   cx0[13] = 61.5; cy0[13] = 55.; cz0[13] = -766.;
   tx0[13] = -44.; ty0[13] = 8.; tz0[13] = -82.;
   //t3_209
   cx1[13] = 77.; cy1[13] = 110.; cz1[13] = -738.;
   tx1[13] = -43.5; ty1[13] = 13.; tz1[13] = -74.;
   //t3_223
   cx15[13] = 682.; cy15[13] = 62.; cz15[13] = -729.;
   tx15[13] = -135.; ty15[13] = -8.; tz15[13] = 99.;

   //t3_224
   cx0[14] = 61.5; cy0[14] = 55.; cz0[14] = -822.;
   tx0[14] = -44.; ty0[14] = 8.; tz0[14] = -82.;
   //t3_225
   cx1[14] = 77.; cy1[14] = 110.; cz1[14] = -792.;
   tx1[14] = -43.5; ty1[14] = 13.; tz1[14] = -74.;
   //t3_239
   cx15[14] = 682.; cy15[14] = 62.; cz15[14] = -783.;
   tx15[14] = -135.; ty15[14] = -8.; tz15[14] = 99.;

   //t3_240
   cx0[15] = 61.5; cy0[15] = 55.; cz0[15] = -879.;
   tx0[15] = -45.; ty0[15] = 8.; tz0[15] = -82.;
   //t3_241
   cx1[15] = 77.; cy1[15] = 110.; cz1[15] = -848.;
   tx1[15] = -43.5; ty1[15] = 13.; tz1[15] = -74.;
   //t3_255
   cx15[15] = 682.; cy15[15] = 62.; cz15[15] = -837.;
   tx15[15] = -135.; ty15[15] = -8.; tz15[15] = 99.;

   tree2->Fill();
   tree2->Write();
   tree3->Fill();
   tree3->Write();

   TVector3 v;
   TGraph2D *pm3[16];
   Int_t nm3;

   Long64_t Ng = 0;
   /*for(Int_t i=0; i<16; i++){
      nm3 = 0;
      f->cd();
      t0[i]->SetBranchAddress("x",&x);
      t0[i]->SetBranchAddress("y",&y);
      t0[i]->SetBranchAddress("z",&z);
      n = t0[i]->GetEntries();
      pm3[i] = new TGraph2D();
      for(Int_t j=0; j<n; j++){
         //if(j%3){continue;}
         t0[i]->GetEntry(j);
         x = x - cx0[i];
         y = y - cy0[i];
         z = z - cz0[i];
         v = TVector3(x,y,z);
         v.RotateX(tx0[i]*DegToRad());
         v.RotateY(ty0[i]*DegToRad());
         v.RotateZ(tz0[i]*DegToRad());
         x = v.X(); y = v.Y(); z = v.Z();
         pm3[i]->SetPoint(nm3++,x,y,z);
      }
   }*/

   ///////////////////////////////////////////////////////////////////////////
   Int_t nr1 = 1000;
   Int_t nr2 = 3000;
   Int_t nr3 = 10000;
   
   TRandom *R = new TRandom();
   R->SetSeed();

   const Double_t E = 0.025 * 2;
   TVector3 ve;
   Long64_t Ng5 = 0;

   TPolyMarker3D *gr0 = new TPolyMarker3D();
   TPolyMarker3D *gr5 = new TPolyMarker3D();

   //側面・小(x<0)
   for(Int_t i=0; i<nr1; i++){
      x = -hX;
      y = -hY + 2*hY*(R->Rndm());
      z = -hZ + 2*hZ*(R->Rndm());
      gr0->SetPoint(Ng++,x,y,z);
   }
   //上面(USならy<0,DSだとy>0)
   for(Int_t i=0; i<nr2; i++){
      y = -hY;
      x = -hX + 2*hX*(R->Rndm());
      z = -hZ + 2*hZ*(R->Rndm());
      gr0->SetPoint(Ng++,x,y,z);
   }
   //側面・大(z>0)
   for(Int_t i=0; i<nr3; i++){
      y = -hY + 2*hY*(R->Rndm());
      x = -hX + 2*hX*(R->Rndm());
      z = hZ;
      gr0->SetPoint(Ng++,x,y,z);
   }
   gr0->SetMarkerStyle(20);
   gr0->SetMarkerSize(0.15);
   gr0->SetMarkerColor(kMagenta);
   gr0->SetName("gr0");
//////////////////////////////////
   //側面・小(x<0)
   for(Int_t i=0; i<nr1; i++){
      x = -hX;
      y = -hY + 2*hY*(R->Rndm());
      z = -hZ + 2*hZ*(R->Rndm());
      y *= 5/4;
      gr5->SetPoint(Ng5++,x,y,z);
   }
   //上面(USならy<0,DSだとy>0)
   for(Int_t i=0; i<nr2; i++){
      y = -hY;
      x = -hX + 2*hX*(R->Rndm());
      z = -hZ + 2*hZ*(R->Rndm());
      y *= 5/4;
      gr5->SetPoint(Ng5++,x,y,z);
   }
   //側面・大(z>0)
   for(Int_t i=0; i<nr3; i++){
      y = -hY + 2*hY*(R->Rndm());
      x = -hX + 2*hX*(R->Rndm());
      z = hZ;
      y *= 5/4;
      gr5->SetPoint(Ng5++,x,y,z);
   }


   gr5->SetMarkerStyle(20);
   gr5->SetMarkerSize(0.15);
   gr5->SetMarkerColor(kCyan);
   gr5->SetName("gr5");
   /*TFile *fb = new TFile("boxmake1.root");
   TTree *tr0 = (TTree*)fb->Get("tr0");
   tr0->SetBranchAddress("x",&x);
   tr0->SetBranchAddress("y",&y);
   tr0->SetBranchAddress("z",&z);
   Int_t Nr0 = tr0->GetEntries();
   TGraph2D *g2d = new TGraph2D();
   for(Int_t i=0; i<Nr0; i++){
      if(i%3){continue;}
      tr0->GetEntry(i);
      g2d->SetPoint(Ng++,x,y,z);
   }
   TCanvas *c[16];
   for(Int_t i=8; i<10; i++){
      c[i] = new TCanvas(Form("counter_%d",i),Form("counter_%d",i),800,600);
      pm3[i]->SetTitle(Form("counter_%d_prefitting; x axis; y axis; z axis",256+i*16));
      pm3[i]->SetMarkerStyle(20);
      pm3[i]->SetMarkerSize(0.15);
      pm3[i]->SetMarkerColor(kBlack);
      pm3[i]->Draw();
      if(i >= 2 && i <= 7){gr5->Draw("same");}
      else{gr0->Draw("same");}
   }*/

   fo->cd();
   gr0->Write();
   gr5->Write();


   TGraph2D *gg = new TGraph2D();
   Ng = 0;
   Int_t k=14;
   Int_t m=15;
   if(m==0){
      t0[k]->SetBranchAddress("x",&x);
      t0[k]->SetBranchAddress("y",&y);
      t0[k]->SetBranchAddress("z",&z);
      gg->SetNameTitle("gg",Form("(%5.2lf,%5.2lf,%5.2lf)cen_(%5.2lf,%5.2lf,%5.2lf)deg; x; y; z",cx0[k],cy0[k],cz0[k],tx0[k],ty0[k],tz0[k]));
      n = t0[k]->GetEntries();
      for(Int_t i=0; i<n; i++){
         if(i%10){continue;}
         t0[k]->GetEntry(i);
         x = x - cx0[k];
         y = y - cy0[k];
         z = z - cz0[k];
         v = TVector3(x,y,z);
         v.RotateX(tx0[k]*DegToRad());
         v.RotateY(ty0[k]*DegToRad());
         v.RotateZ(tz0[k]*DegToRad());
         x = v.X(); y = v.Y(); z = v.Z();
         gg->SetPoint(Ng++,x,y,z);
      }
   }
   else if(m==1){
      t1[k]->SetBranchAddress("x",&x);
      t1[k]->SetBranchAddress("y",&y);
      t1[k]->SetBranchAddress("z",&z);
      gg->SetNameTitle("gg",Form("(%5.2lf,%5.2lf,%5.2lf)cen_(%5.2lf,%5.2lf,%5.2lf)deg; x; y; z",cx1[k],cy1[k],cz1[k],tx1[k],ty1[k],tz1[k]));
      n = t1[k]->GetEntries();
      for(Int_t i=0; i<n; i++){
         if(i%10){continue;}
         t1[k]->GetEntry(i);
         x = x - cx1[k];
         y = y - cy1[k];
         z = z - cz1[k];
         v = TVector3(x,y,z);
         v.RotateX(tx1[k]*DegToRad());
         v.RotateY(ty1[k]*DegToRad());
         v.RotateZ(tz1[k]*DegToRad());
         x = v.X(); y = v.Y(); z = v.Z();
         gg->SetPoint(Ng++,x,y,z);
      }
   }
   else if(m==15){
      t15[k]->SetBranchAddress("x",&x);
      t15[k]->SetBranchAddress("y",&y);
      t15[k]->SetBranchAddress("z",&z);
      gg->SetNameTitle("gg",Form("(%5.2lf,%5.2lf,%5.2lf)cen_(%5.2lf,%5.2lf,%5.2lf)deg; x; y; z",
               cx15[k],cy15[k],cz15[k],tx15[k],ty15[k],tz15[k]));
      n = t15[k]->GetEntries();
      for(Int_t i=0; i<n; i++){
         if(i%10){continue;}
         t15[k]->GetEntry(i);
         x = x - cx15[k];
         y = y - cy15[k];
         z = z - cz15[k];
         v = TVector3(x,y,z);
         v.RotateX(tx15[k]*DegToRad());
         v.RotateY(ty15[k]*DegToRad());
         v.RotateZ(tz15[k]*DegToRad());
         x = v.X(); y = v.Y(); z = v.Z();
         gg->SetPoint(Ng++,x,y,z);
      }
   }

   TCanvas *c1 = new TCanvas("c1","c1",0,0,1200,400);
   c1->Divide(2,1);
   c1->cd(1);
   if(m==0){t0[k]->Draw("x:y:z");}
   else if(m==1){t1[k]->Draw("x:y:z");}
   else if(m==15){t15[k]->Draw("x:y:z");}
   c1->cd(2);
   gg->SetMarkerStyle(20);
   gg->SetMarkerSize(0.4);
   gg->Draw();
}

