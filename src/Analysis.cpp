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

void Analysis()
{
   Double_t cx,cy,cz,tx,ty,tz;
   Double_t cx2,cy2,cz2,tx2,ty2,tz2;
   Double_t ncos;
   TVector3 cv0[16],n0[16],cv1[16],n1[16];
   TVector3 cv15[16],n15[16];
   TVector3 top0[16],top1[16],top15[16];
   Double_t dcv[45],dtop0_15[16],dtop0_1[16],dtop1_15[16],dncos[45];
   TVector3 nCAD0[16],nCAD1[16],nCAD15[16];
   TVector3 cmCAD0[16],cmCAD1[16],cmCAD15[16];
   TVector3 topCAD0[16],topCAD1[16],topCAD15[16];
   TVector3 v,v_adp,cvd,top_adp;
   Double_t dcvCAD[45],dncosCAD[45],dtopCAD0_15[16],dtopCAD0_1[16],dtopCAD1_15[16];
   TVector3 cHole0[23],cHole1[18],cHole15[21];
   Double_t dcvHole[45],dcvHole0_15[16];
   TFile *f,*f1,*f15;
   TTree *t,*t1,*t15;
   Int_t j,k;
   for(Int_t i=0; i<16; i++){
      j = i * 16 + 256;
      f = new TFile(Form("data/dir_migrad/migrad_%d.root",j));
      t = (TTree*)f->Get("par2tree");
      t->SetBranchAddress("cx",&cx); t->SetBranchAddress("cy",&cy);
      t->SetBranchAddress("cz",&cz); t->SetBranchAddress("tx",&tx);
      t->SetBranchAddress("ty",&ty); t->SetBranchAddress("tz",&tz);
      t->SetBranchAddress("cx2",&cx2); t->SetBranchAddress("cy2",&cy2);
      t->SetBranchAddress("cz2",&cz2); t->SetBranchAddress("tx2",&tx2);
      t->SetBranchAddress("ty2",&ty2); t->SetBranchAddress("tz2",&tz2);
      t->GetEntry(0);

      if(i>=8 && i<=13){
         v_adp = TVector3(-1.65,5.,0.5); //カバーの厚み
         v_adp.RotateZ(-tz2 * DegToRad()); v_adp.RotateY(-ty2 * DegToRad());
         v_adp.RotateX(-tx2 * DegToRad()); v_adp.RotateZ(-tz * DegToRad());
         v_adp.RotateY(-ty * DegToRad()); v_adp.RotateX(-tx * DegToRad());
      }else{
         v_adp = TVector3(-1.65,0.,0.5);
         v_adp.RotateZ(-tz2 * DegToRad()); v_adp.RotateY(-ty2 * DegToRad());
         v_adp.RotateX(-tx2 * DegToRad()); v_adp.RotateZ(-tz * DegToRad());
         v_adp.RotateY(-ty * DegToRad()); v_adp.RotateX(-tx * DegToRad());
      }

      v = TVector3(cx2,cy2,cz2);
      v.RotateZ(-tz * DegToRad()); v.RotateY(-ty * DegToRad()); v.RotateX(-tx * DegToRad());
      cv0[i] = TVector3(cx,cy,cz);

      top_adp = TVector3(-1.65,-20.,0.5);
      top_adp.RotateZ(-tz2 * DegToRad()); top_adp.RotateY(-ty2 * DegToRad());
      top_adp.RotateX(-tx2 * DegToRad()); top_adp.RotateZ(-tz * DegToRad());
      top_adp.RotateY(-ty * DegToRad()); top_adp.RotateX(-tx * DegToRad());
      top0[i] = cv0[i] + v + top_adp;
      cv0[i] = cv0[i] + v + v_adp;

      n0[i] = TVector3(0,0,1);
      n0[i].RotateZ(-tz2 * DegToRad()); n0[i].RotateY(-ty2 * DegToRad());
      n0[i].RotateX(-tx2 * DegToRad()); n0[i].RotateZ(-tz * DegToRad());
      n0[i].RotateY(-ty * DegToRad()); n0[i].RotateX(-tx * DegToRad());
      if(i >= 1 ){ 
         cvd = cv0[i-1] - cv0[i];
         dcv[i-1] = cvd.Mag();
         dncos[i-1] = ACos(n0[i]*n0[0])*RadToDeg();
      }else{continue;}
   }
   TGraph2D *gdeg0 = new TGraph2D();
   gdeg0->SetNameTitle("gdeg0",
         "Angle b/w normal vectors BP0 from fitting; No.; No.; Angle [Degree]");
   for(Int_t i=0; i<16; i++){
      for(Int_t m=0; m<16; m++){
         k = 16 * i + m;
         ncos = ACos(n0[i]*n0[m])*RadToDeg();
         gdeg0->SetPoint(k,i,m,ncos);
      }
   }
   for(Int_t i=0; i<16; i++){
      j = i * 16 + 257;
      f1 = new TFile(Form("data/dir_migrad/migrad_%d.root",j));
      t1 = (TTree*)f1->Get("par2tree");
      t1->SetBranchAddress("cx",&cx); t1->SetBranchAddress("cy",&cy);
      t1->SetBranchAddress("cz",&cz); t1->SetBranchAddress("tx",&tx);
      t1->SetBranchAddress("ty",&ty); t1->SetBranchAddress("tz",&tz);
      t1->SetBranchAddress("cx2",&cx2); t1->SetBranchAddress("cy2",&cy2);
      t1->SetBranchAddress("cz2",&cz2); t1->SetBranchAddress("tx2",&tx2);
      t1->SetBranchAddress("ty2",&ty2); t1->SetBranchAddress("tz2",&tz2);
      t1->GetEntry(0);

      if(i>=7 && i<=12){
         v_adp = TVector3(-1.65,5.,0.5);
         v_adp.RotateZ(-tz2 * DegToRad()); v_adp.RotateY(-ty2 * DegToRad());
         v_adp.RotateX(-tx2 * DegToRad()); v_adp.RotateZ(-tz * DegToRad());
         v_adp.RotateY(-ty * DegToRad()); v_adp.RotateX(-tx * DegToRad());
      }else{
         v_adp = TVector3(-1.65,0.,0.5);
         v_adp.RotateZ(-tz2 * DegToRad()); v_adp.RotateY(-ty2 * DegToRad());
         v_adp.RotateX(-tx2 * DegToRad()); v_adp.RotateZ(-tz * DegToRad());
         v_adp.RotateY(-ty * DegToRad()); v_adp.RotateX(-tx * DegToRad());
      }

      v = TVector3(cx2,cy2,cz2);
      v.RotateZ(-tz * DegToRad());
      v.RotateY(-ty * DegToRad());
      v.RotateX(-tx * DegToRad());
      cv1[i] = TVector3(cx,cy,cz);
      top_adp = TVector3(-1.65,-20.,0.5);
      top_adp.RotateZ(-tz2 * DegToRad()); top_adp.RotateY(-ty2 * DegToRad());
      top_adp.RotateX(-tx2 * DegToRad()); top_adp.RotateZ(-tz * DegToRad());
      top_adp.RotateY(-ty * DegToRad()); top_adp.RotateX(-tx * DegToRad());
      top1[i] = cv1[i] + v + top_adp;
      cvd = top1[i] - top0[i];
      dtop0_1[i] = cvd.Mag();
      cv1[i] = cv1[i] + v + v_adp;

      n1[i] = TVector3(0,0,1);
      n1[i].RotateZ(-tz2 * DegToRad());
      n1[i].RotateY(-ty2 * DegToRad());
      n1[i].RotateX(-tx2 * DegToRad());
      n1[i].RotateZ(-tz * DegToRad());
      n1[i].RotateY(-ty * DegToRad());
      n1[i].RotateX(-tx * DegToRad());
      if(i >= 1 ){ 
         cvd = cv1[i-1] - cv1[i];
         dcv[i+14] = cvd.Mag();
         dncos[i+14] = ACos(n1[i]*n1[0])*RadToDeg();
      }else{continue;}
   }
   TGraph2D *gdeg1 = new TGraph2D();
   gdeg1->SetNameTitle("gdeg1",
         "Angle b/w normal vectors BP1 from fitting; No.; No.; Angle [Degree]");
   for(Int_t i=0; i<16; i++){
      for(Int_t m=0; m<16; m++){
         k = 16 * i + m;
         ncos = ACos(n1[i]*n1[m])*RadToDeg();
         gdeg1->SetPoint(k,i,m,ncos);
      }
   }
   for(Int_t i=0; i<16; i++){
      j = i * 16 + 271;
      f15 = new TFile(Form("data/dir_migrad/migrad_%d.root",j));
      t15 = (TTree*)f15->Get("par2tree");
      t15->SetBranchAddress("cx",&cx); t15->SetBranchAddress("cy",&cy);
      t15->SetBranchAddress("cz",&cz); t15->SetBranchAddress("tx",&tx);
      t15->SetBranchAddress("ty",&ty); t15->SetBranchAddress("tz",&tz);
      t15->SetBranchAddress("cx2",&cx2); t15->SetBranchAddress("cy2",&cy2);
      t15->SetBranchAddress("cz2",&cz2); t15->SetBranchAddress("tx2",&tx2);
      t15->SetBranchAddress("ty2",&ty2); t15->SetBranchAddress("tz2",&tz2);
      t15->GetEntry(0);

      if(i>=7 && i<=12){
         v_adp = TVector3(-1.65,5.,0.5);
         v_adp.RotateZ(-tz2 * DegToRad()); v_adp.RotateY(-ty2 * DegToRad());
         v_adp.RotateX(-tx2 * DegToRad()); v_adp.RotateZ(-tz * DegToRad());
         v_adp.RotateY(-ty * DegToRad()); v_adp.RotateX(-tx * DegToRad());
      }else{
         v_adp = TVector3(-1.65,0.,0.5);
         v_adp.RotateZ(-tz2 * DegToRad()); v_adp.RotateY(-ty2 * DegToRad());
         v_adp.RotateX(-tx2 * DegToRad()); v_adp.RotateZ(-tz * DegToRad());
         v_adp.RotateY(-ty * DegToRad()); v_adp.RotateX(-tx * DegToRad());
      }

      v = TVector3(cx2,cy2,cz2);
      v.RotateZ(-tz * DegToRad()); v.RotateY(-ty * DegToRad()); v.RotateX(-tx * DegToRad());
      cv15[i] = TVector3(cx,cy,cz);

      top_adp = TVector3(-1.65,-20.,0.5);
      top_adp.RotateZ(-tz2 * DegToRad()); top_adp.RotateY(-ty2 * DegToRad());
      top_adp.RotateX(-tx2 * DegToRad()); top_adp.RotateZ(-tz * DegToRad());
      top_adp.RotateY(-ty * DegToRad()); top_adp.RotateX(-tx * DegToRad());
      top15[i] = cv15[i] + v + top_adp;
      cvd = top15[i] - top0[i];
      dtop0_15[i] = cvd.Mag();
      cvd = top15[i] - top1[i];
      dtop1_15[i] = cvd.Mag();
      cv15[i] = cv15[i] + v + v_adp;

      n15[i] = TVector3(0,0,1);
      n15[i].RotateZ(-tz2 * DegToRad()); n15[i].RotateY(-ty2 * DegToRad());
      n15[i].RotateX(-tx2 * DegToRad()); n15[i].RotateZ(-tz * DegToRad());
      n15[i].RotateY(-ty * DegToRad()); n15[i].RotateX(-tx * DegToRad());

      if(i >= 1 ){ 
         cvd = cv15[i-1] - cv15[i];
         dcv[i+29] = cvd.Mag();
         dncos[i+29] = ACos(n15[i]*n15[0])*RadToDeg();
      }else{continue;}
   }
   TGraph2D *gdeg15 = new TGraph2D();
   gdeg15->SetNameTitle("gdeg15",
         "Angle b/w normal vectors BP15 from fitting; No.; No.; Angle [Degree]");
   for(Int_t i=0; i<16; i++){
      for(Int_t m=0; m<16; m++){
         k = 16 * i + m;
         ncos = ACos(n15[i]*n15[m])*RadToDeg();
         gdeg15->SetPoint(k,i,m,ncos);
      }
   }

   nCAD0[15] = TVector3(-0.185,-3.53075,-3.5); cmCAD0[15] = TVector3(316.575,-16.387375,-1125);
   nCAD0[14] = TVector3(-0.185,-3.53075,-3.5); cmCAD0[14] = TVector3(316.575,-16.387375,-1070);
   nCAD0[13] = TVector3(-0.185,-3.53075,-3.515); cmCAD0[13] = TVector3(316.575,-16.387375,-1015.0125);
   nCAD0[12] = TVector3(-0.185,-3.53075,-3.515); cmCAD0[12] = TVector3(316.575,-16.387375,-960.0125);
   nCAD0[11] = TVector3(-0.185,-3.53075,-3.53); cmCAD0[11] = TVector3(316.575,-16.387375,-905);
   nCAD0[10] = TVector3(-0.185,-3.53075,-3.53); cmCAD0[10] = TVector3(316.575,-16.387375,-850);
   nCAD0[9] = TVector3(-0.185,-3.53075,-3.53); cmCAD0[9] = TVector3(316.575,-16.387375,-795);
   nCAD0[8] = TVector3(-0.185,-3.53075,-3.53); cmCAD0[8] = TVector3(316.575,-16.387375,-740);
   nCAD0[7] = TVector3(-0.1825,-3.53075,-3.53); cmCAD0[7] = TVector3(321.56875,-16.645875,-685);
   nCAD0[6] = TVector3(-0.1825,-3.53075,-3.53); cmCAD0[6] = TVector3(321.56875,-16.645875,-630);
   nCAD0[5] = TVector3(-0.1825,-3.53075,-3.53); cmCAD0[5] = TVector3(321.56875,-16.645875,-575);
   nCAD0[4] = TVector3(-0.1825,-3.53075,-3.53); cmCAD0[4] = TVector3(321.56875,-16.645875,-520);
   nCAD0[3] = TVector3(-0.1825,-3.53075,-3.53); cmCAD0[3] = TVector3(321.56875,-16.645875,-465);
   nCAD0[2] = TVector3(-0.1825,-3.53075,-3.53); cmCAD0[2] = TVector3(321.56875,-16.645875,-410);
   nCAD0[1] = TVector3(-0.185,-3.53075,-3.53); cmCAD0[1] = TVector3(316.575,-16.387375,-355);
   nCAD0[0] = TVector3(-0.185,-3.53075,-3.53); cmCAD0[0] = TVector3(316.575,-16.387375,-300);

   nCAD1[15] = TVector3(-0.8125,-3.443,-3.5);
   cmCAD1[15] = TVector3(308.52625,-72.7995,-1097.5);
   nCAD1[14] = TVector3(-0.8125,-3.443,-3.495);
   cmCAD1[14] = TVector3(308.52625,-72.7995,-1042.5025);
   nCAD1[13] = TVector3(-0.8125,-3.443,-3.515);
   cmCAD1[13] = TVector3(308.52625,-72.7995,-987.5125);
   nCAD1[12] = TVector3(-0.8125,-3.443,-3.53);
   cmCAD1[12] = TVector3(308.52625,-72.7995,-932.5);
   nCAD1[11] = TVector3(-0.8125,-3.443,-3.53);
   cmCAD1[11] = TVector3(308.52625,-72.7995,-877.5);
   nCAD1[10] = TVector3(-0.8125,-3.443,-3.53);
   cmCAD1[10] = TVector3(308.52625,-72.7995,-822.5);
   nCAD1[9] = TVector3(-0.8125,-3.443,-3.53);
   cmCAD1[9] = TVector3(308.52625,-72.7995,-767.5);
   nCAD1[8] = TVector3(-0.81,-3.4405,-3.53);
   cmCAD1[8] = TVector3(313.3925,-73.9475,-712.5);
   nCAD1[7] = TVector3(-0.81,-3.4405,-3.53);
   cmCAD1[7] = TVector3(313.3925,-73.9475,-657.5);
   nCAD1[6] = TVector3(-0.81,-3.4405,-3.53);
   cmCAD1[6] = TVector3(313.3925,-73.9475,-602.5);
   nCAD1[5] = TVector3(-0.81,-3.4405,-3.53);
   cmCAD1[5] = TVector3(313.3925,-73.9475,-547.5);
   nCAD1[4] = TVector3(-0.81,-3.4405,-3.53);
   cmCAD1[4] = TVector3(313.3925,-73.9475,-492.5);
   nCAD1[3] = TVector3(-0.81,-3.4405,-3.53);
   cmCAD1[3] = TVector3(313.3925,-73.9475,-437.5);
   nCAD1[2] = TVector3(-0.8125,-3.443,-3.53);
   cmCAD1[2] = TVector3(308.52625,-72.7995,-382.5);
   nCAD1[1] = TVector3(-0.8125,-3.443,-3.53);
   cmCAD1[1] = TVector3(308.52625,-72.7995,-327.5);
   nCAD1[0] = TVector3(-0.8125,-3.443,-3.53);
   cmCAD1[0] = TVector3(308.52625,-72.7995,-272.5);
   
   nCAD15[15] = TVector3(-1.345,3.27,-3.5);
   cmCAD15[15] = TVector3(-293.2125,-120.48225,-1097.5);
   nCAD15[14] = TVector3(-1.345,3.27,-3.495);
   cmCAD15[14] = TVector3(-293.2125,-120.48225,-1042.5025);
   nCAD15[13] = TVector3(-1.345,3.27,-3.515);
   cmCAD15[13] = TVector3(-293.2125,-120.48225,-987.5125);
   nCAD15[12] = TVector3(-1.345,3.27,-3.53);
   cmCAD15[12] = TVector3(-293.2125,-120.48225,-932.5);
   nCAD15[11] = TVector3(-1.345,3.27,-3.53);
   cmCAD15[11] = TVector3(-293.2125,-120.48225,-877.5);
   nCAD15[10] = TVector3(-1.345,3.27,-3.53);
   cmCAD15[10] = TVector3(-293.2125,-120.48225,-822.5);
   nCAD15[9] = TVector3(-1.345,3.27,-3.53);
   cmCAD15[9] = TVector3(-293.2125,-120.48225,-767.5);
   nCAD15[8] = TVector3(-1.345,3.27025,-3.53);
   cmCAD15[8] = TVector3(-297.8375,-122.382375,-712.5);
   nCAD15[7] = TVector3(-1.345,3.27025,-3.53);
   cmCAD15[7] = TVector3(-297.8375,-122.382375,-657.5);
   nCAD15[6] = TVector3(-1.345,3.27025,-3.53);
   cmCAD15[6] = TVector3(-297.8375,-122.382375,-602.5);
   nCAD15[5] = TVector3(-1.345,3.27025,-3.53);
   cmCAD15[5] = TVector3(-297.8375,-122.382375,-547.5);
   nCAD15[4] = TVector3(-1.345,3.27025,-3.53);
   cmCAD15[4] = TVector3(-297.8375,-122.382375,-492.5);
   nCAD15[3] = TVector3(-1.345,3.27025,-3.53);
   cmCAD15[3] = TVector3(-297.8375,-122.382375,-437.5);
   nCAD15[2] = TVector3(-1.345,3.27,-3.53);
   cmCAD15[2] = TVector3(-293.2125,-120.48225,-382.5);
   nCAD15[1] = TVector3(-1.345,3.27,-3.53);
   cmCAD15[1] = TVector3(-293.2125,-120.48225,-327.5);
   nCAD15[0] = TVector3(-1.345,3.27,-3.53);
   cmCAD15[0] = TVector3(-293.2125,-120.48225,-272.5);

   topCAD0[15] = TVector3(296.6025,-15.3535,-1125);
   topCAD0[14] = TVector3(296.6025,-15.3535,-1070);
   topCAD0[13] = TVector3(296.6025,-15.3535,-1015.0125);
   topCAD0[12] = TVector3(296.6025,-15.3535,-960.0125);
   topCAD0[11] = TVector3(296.6025,-15.3535,-905);
   topCAD0[10] = TVector3(296.6025,-15.3535,-850);
   topCAD0[9] = TVector3(296.6025,-15.3535,-795);
   topCAD0[8] = TVector3(296.6025,-15.3535,-740);
   topCAD0[7] = TVector3(296.6025,-15.3535,-685);
   topCAD0[6] = TVector3(296.6025,-15.3535,-630);
   topCAD0[5] = TVector3(296.6025,-15.3535,-575);
   topCAD0[4] = TVector3(296.6025,-15.3535,-520);
   topCAD0[3] = TVector3(296.6025,-15.3535,-465);
   topCAD0[2] = TVector3(296.6025,-15.3535,-410);
   topCAD0[1] = TVector3(296.6025,-15.3535,-355);
   topCAD0[0] = TVector3(296.6025,-15.3535,-300);
   topCAD15[15] = TVector3(-274.715,-112.8815,-1097.5);
   topCAD15[14] = TVector3(-274.715,-112.8815,-1042.5025);
   topCAD15[13] = TVector3(-274.715,-112.8815,-987.5125);
   topCAD15[12] = TVector3(-274.715,-112.8815,-932.5);
   topCAD15[11] = TVector3(-274.715,-112.8815,-877.5);
   topCAD15[10] = TVector3(-274.715,-112.8815,-822.5);
   topCAD15[9] = TVector3(-274.715,-112.8815,-767.5);
   topCAD15[8] = TVector3(-274.715,-112.8815,-712.5);
   topCAD15[7] = TVector3(-274.715,-112.8815,-657.5);
   topCAD15[6] = TVector3(-274.715,-112.8815,-602.5);
   topCAD15[5] = TVector3(-274.715,-112.8815,-547.5);
   topCAD15[4] = TVector3(-274.715,-112.8815,-492.5);
   topCAD15[3] = TVector3(-274.715,-112.8815,-437.5);
   topCAD15[2] = TVector3(-274.715,-112.8815,-382.5);
   topCAD15[1] = TVector3(-274.715,-112.8815,-327.5);
   topCAD15[0] = TVector3(-274.715,-112.8815,-272.5);
   
   TGraph2D *gdegCAD0 = new TGraph2D();
   gdegCAD0->SetNameTitle("gdegCAD0",
         "Angle b/w normal vectors BP0 from CAD; No.; No.; Angle [Degree]");
   TGraph2D *gdegCAD1 = new TGraph2D();
   gdegCAD1->SetNameTitle("gdegCAD1",
         "Angle b/w normal vectors BP1 from CAD; No.; No.; Angle [Degree]");
   TGraph2D *gdegCAD15 = new TGraph2D();
   gdegCAD15->SetNameTitle("gdegCAD15",
         "Angle b/w normal vectors BP15 from CAD; No.; No.; Angle [Degree]");
   for(Int_t i=0; i<16; i++){
      cvd = topCAD15[i] - topCAD0[i];
      dtopCAD0_15[i] = cvd.Mag();
      nCAD0[i] = nCAD0[i].Unit();
      nCAD1[i] = nCAD1[i].Unit();
      nCAD15[i] = nCAD15[i].Unit();
      if(i >= 1){
         cvd = cmCAD0[i-1] - cmCAD0[i];
         dcvCAD[i-1] = cvd.Mag();
         cvd = cmCAD1[i-1] - cmCAD1[i];
         dcvCAD[i+14] = cvd.Mag();
         cvd = cmCAD15[i-1] - cmCAD15[i];
         dcvCAD[i+29] = cvd.Mag();
         dncosCAD[i-1] = ACos(nCAD0[i]*nCAD0[0])*RadToDeg();
         dncosCAD[i+14] = ACos(nCAD1[i]*nCAD1[0])*RadToDeg();
         dncosCAD[i+29] = ACos(nCAD15[i]*nCAD15[0])*RadToDeg();
      }
   }
   for(Int_t i=0; i<16; i++){
      for(Int_t m=0; m<16; m++){
         k = 16 * i + m;
         ncos = ACos(nCAD0[i]*nCAD0[m])*RadToDeg();
         gdegCAD0->SetPoint(k,i,m,ncos);
         ncos = ACos(nCAD1[i]*nCAD1[m])*RadToDeg();
         gdegCAD1->SetPoint(k,i,m,ncos);
         ncos = ACos(nCAD15[i]*nCAD15[m])*RadToDeg();
         gdegCAD15->SetPoint(k,i,m,ncos);
      }
   }

   cHole0[19] = TVector3(-378.698,-66.0343,-27.0686); //496 center
   cHole0[18] = TVector3(-378.655,-66.0343,27.87635); //480 center
   cHole0[17] = TVector3(-378.615,-66.01705,83.02945); //464 center
   cHole0[16] = TVector3(-378.591,-66.00735,137.9032); //448 center
   cHole0[15] = TVector3(-378.6375,-66.08215,192.8965); //432 center
   cHole0[14] = TVector3(-378.6615,-66.0271,247.8325); //416 center
   cHole0[13] = TVector3(-378.766,-66.05105,302.8995); //400 center
   cHole0[12] = TVector3(-378.622,-66.05465,357.8455); //384 center
   ////
   cHole0[11] = TVector3(-385.523,-116.047,307.32); //384 ch2
   cHole0[10] = TVector3(-385.581,-116.256,362.352); //368 ch2
   cHole0[9] = TVector3(-385.796,-116.156,417.279); //352 ch2
   ////
   cHole0[8] = TVector3(-378.781,-66.0805,467.8615); //352 center
   cHole0[7] = TVector3(-378.621,-66.0014,522.8505); //336 center

   cHole0[19] = TVector3(-378.698,-66.0343,-27.0686); //496 center
   cHole0[18] = TVector3(-378.655,-66.0343,27.87635); //480 center
   cHole0[17] = TVector3(-378.615,-66.01705,83.02945); //464 center
   cHole0[16] = TVector3(-378.591,-66.00735,137.9032); //448 center
   cHole0[15] = TVector3(-378.6375,-66.08215,192.8965); //432 center
   cHole0[14] = TVector3(-378.6615,-66.0271,247.8325); //416 center
   cHole0[13] = TVector3(-378.766,-66.05105,302.8995); //400 center
   cHole0[12] = TVector3(-378.622,-66.05465,357.8455); //384 center
   ////
   cHole0[11] = TVector3(-385.523,-116.047,307.32); //384 ch2
   cHole0[10] = TVector3(-385.581,-116.256,362.352); //368 ch2
   cHole0[9] = TVector3(-385.796,-116.156,417.279); //352 ch2
   ////
   cHole0[8] = TVector3(-378.781,-66.0805,467.8615); //352 center
   cHole0[7] = TVector3(-378.621,-66.0014,522.8505); //336 center
   cHole0[6] = TVector3(-379.0855,-65.9106,577.7875); //320 center
   cHole0[5] = TVector3(-378.87,-66.02615,632.798); //304 center
   cHole0[4] = TVector3(-378.6035,-66.2147,687.959); //288 center
   ////
   cHole0[3] = TVector3(-385.725,-116.239,637.287); //288 ch2
   cHole0[2] = TVector3(-385.14,-116.168,692.708);  //272 ch2
   ////
   cHole0[1] = TVector3(-337.161,-21.5341,792.153);  //272 ch1
   cHole0[0] = TVector3(-337.378,-21.6332,847.085); //256 ch1

   cHole0[20] = TVector3(-371.482,-16.1904,738.631); //288 ch1
   cHole0[21] = TVector3(-371.783,-15.8288,573.412); //336 ch1
   cHole0[22] = TVector3(-372.8145,-70.58865,409.1815); //368 center
   for(Int_t i=0; i<15; i++){
      if(i==0 || i==1) dcvHole[i] = (cHole0[2*i+1]-cHole0[2*i]).Mag();
      else if(i>=2 && i<=5) dcvHole[i] = (cHole0[i+3]-cHole0[i+2]).Mag();
      else if(i>=6 && i<=7) dcvHole[i] = (cHole0[i+4]-cHole0[i+3]).Mag();
      else if(i>=8 && i<=14) dcvHole[i] = (cHole0[i+5]-cHole0[i+4]).Mag();
   }


   cHole1[17] = TVector3(-363.2995,1.10275,0.5409); //497center
   cHole1[16] = TVector3(-363.476,0.8724,55.39474); //481center
   cHole1[15] = TVector3(-363.495,0.78445,110.4411); //465center
   cHole1[14] = TVector3(-363.1625,1.07355,165.423); //449center
   cHole1[13] = TVector3(-363.2425,0.9732,220.4715); //433center
   cHole1[12] = TVector3(-363.346,0.97715,275.3945); //417center
   cHole1[11] = TVector3(-363.1605,0.91005,330.418); //401center
   cHole1[10] = TVector3(-363.5115,0.94785,385.367); //385center
   cHole1[9] = TVector3(-363.359,0.8175,440.358); //369center
   cHole1[8] = TVector3(-363.526,0.98945,495.355); //353center
   cHole1[7] = TVector3(-363.494,0.8519,550.286); //337center
   cHole1[6] = TVector3(-363.3485,0.937,605.3795); //321center
   cHole1[5] = TVector3(-363.4925,0.8751,660.365); //305center
   ////
   cHole1[4] = TVector3(-379.278,-47.1016,609.815); //305ch2
   cHole1[3] = TVector3(-379.043,-47.1232,664.826); //289ch2
   cHole1[2] = TVector3(-379.459,-47.2953,719.955); //273ch2
   ////
   cHole1[1] = TVector3(-314.339,37.3036,819.668); //273ch1
   cHole1[0] = TVector3(-314.437,37.2272,874.588); //257ch1
   for(Int_t i=0; i<15; i++){
      if(i==0) dcvHole[i+15] = (cHole1[i+1]-cHole1[i]).Mag();
      else if(i>=1 && i<=2) dcvHole[i+15] = (cHole1[i+2]-cHole1[i+1]).Mag();
      else if(i>=3 && i<=14) dcvHole[i+15] = (cHole1[i+3]-cHole1[i+2]).Mag();
   }

   cHole15[19] = TVector3(361.963,0.9977,0.34075); //511center
   cHole15[18] = TVector3(362.0175,0.95105,55.26055); //495center
   cHole15[17] = TVector3(361.9255,0.92905,110.50445); //479center
   cHole15[16] = TVector3(362.0215,0.85475,165.5565); //463center
   cHole15[15] = TVector3(362.0765,0.71475,220.545); //447center
   cHole15[14] = TVector3(362.142,0.7479,275.482); //431center
   cHole15[13] = TVector3(362.245,0.6988,330.489); //415center
   cHole15[12] = TVector3(362.162,0.67075,385.4505); //399center
   cHole15[11] = TVector3(362.406,0.68585,440.37); //383center
   cHole15[10] = TVector3(362.5515,0.65715,495.4065); //367center
   ////
   cHole15[9] = TVector3(378.42,-47.4245,545.486); //367ch1
   cHole15[8] = TVector3(378.506,-47.4444,600.389); //351ch1
   cHole15[7] = TVector3(378.799,-47.5651,655.46); //335ch1
   ////
   cHole15[6] = TVector3(362.76,0.6287,605.4165); //335center
   cHole15[5] = TVector3(362.998,0.7015,660.2235); //319center
   ////
   cHole15[4] = TVector3(347.401,48.8483,610.045); //319ch2
   cHole15[3] = TVector3(347.382,48.7188,665.039); //303ch2
   cHole15[2] = TVector3(347.359,48.8434,720.112); //287ch2
   ////
   cHole15[1] = TVector3(343.275,-59.6307,819.773); //287ch1
   cHole15[0] = TVector3(343.466,-59.3359,874.658); //271ch1
   cHole15[20] = TVector3(343.82,-59.3295,764.857); //303ch1
   for(Int_t i=0; i<15; i++){
      if(i==0) dcvHole[i+30] = (cHole15[i+1]-cHole15[i]).Mag();
      else if(i>=1 && i<=2) dcvHole[i+30] = (cHole15[i+2]-cHole15[i+1]).Mag();
      else if(i==3) dcvHole[i+30] = (cHole15[i+3]-cHole15[i+2]).Mag();
      else if(i>=4 && i<=5) dcvHole[i+30] = (cHole15[i+4]-cHole15[i+3]).Mag();
      else if(i>=6 && i<=14) dcvHole[i+30] = (cHole15[i+5]-cHole15[i+4]).Mag();
   }
   dcvHole0_15[0] = (cHole15[0]-cHole0[0]).Mag();
   dcvHole0_15[1] = (cHole15[1]-cHole0[1]).Mag();
   dcvHole0_15[2] = (cHole15[20]-cHole0[20]).Mag();
   dcvHole0_15[3] = (cHole15[5]-cHole0[5]).Mag();
   dcvHole0_15[4] = (cHole15[6]-cHole0[6]).Mag();
   dcvHole0_15[5] = (cHole15[8]-cHole0[21]).Mag();
   dcvHole0_15[6] = (cHole15[10]-cHole0[8]).Mag();
   dcvHole0_15[7] = (cHole15[11]-cHole0[22]).Mag();

   for(Int_t i=8; i<16; i++){
      dcvHole0_15[i] = (cHole15[i+4]-cHole0[i+4]).Mag();
   }


   std::ofstream ofs("/Users/taku/Documents/survey/outputfiles/Analysis.csv");
   for(Int_t i=0; i<16; i++){
      ofs << cv0[i].X() << ", " << cv0[i].Y() << ", " << cv0[i].Z() << std::endl;
   }
   for(Int_t i=0; i<16; i++){
      ofs << cv1[i].X() << ", " << cv1[i].Y() << ", " << cv1[i].Z() << std::endl;
   }
   for(Int_t i=0; i<16; i++){
      ofs << cv15[i].X() << ", " << cv15[i].Y() << ", " << cv15[i].Z() << std::endl;
   }
   ofs.close();

   /*TH1D *i1 = new TH1D("i1","Intervals [mm] (Fitting)",100,50,60);
   TH1D *i2 = new TH1D("i2","Intervals [mm] (CAD)",100,50,60);
   TH1D *i3 = new TH1D("i3","Intervals [mm] (Holes)",100,50,60);
   TH1D *e1 = new TH1D("e1","Edge intervals [mm] (Fitting)",100,570,590);
   TH1D *e2 = new TH1D("e2","Edge intervals [mm] (CAD)",100,570,590);
   TH1D *a1 = new TH1D("a1","Angles b/w n&n0 [Degree] (Fitting)",100,-3,3);
   TH1D *a2 = new TH1D("a2","Angles b/w n&n0 [Degree] (CAD)",100,-3,3);
   TFile *fout = new TFile("Analysis.root","recreate");
   TTree *tout = new TTree("tout","tout");
   tout->Branch("dcv",dcv,"dcv[45]/D");
   tout->Branch("dcvCAD",dcvCAD,"dcvCAD[45]/D");
   tout->Branch("dcvHole",dcvHole,"dcvHole[45]/D");
   tout->Branch("Angle",dncos,"Angle[45]/D");
   tout->Branch("AngleCAD",dncosCAD,"AngleCAD[45]/D");
   tout->Branch("dEdge",dtop0_15,"dEdge[16]/D");
   tout->Branch("dEdgeCAD",dtopCAD0_15,"dEdgeCAD[16]/D");
   tout->Branch("dEdgeHole",dcvHole0_15,"dEdgeHole[16]/D");
   tout->Fill();
   tout->Write();
   gdeg0->Write();
   gdeg1->Write();
   gdeg15->Write();
   gdegCAD0->Write();
   gdegCAD1->Write();
   gdegCAD15->Write();
   i1->Write();
   i2->Write();
   i3->Write();
   e1->Write();
   e2->Write();
   a1->Write();
   a2->Write();
   */
}
