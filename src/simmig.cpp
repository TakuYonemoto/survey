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

void Di(Double_t x, Double_t y, Double_t z, Double_t &r);
void Di2(Double_t x, Double_t y, Double_t z, Double_t &r);
void Loss(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void first_fit();
void refit();
void rerefit();
void exfit_yaxis1();
void exfit_yaxis2();
void adjfit_yaxis();
void adjfit_xzaxis();
//void rerefit();

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
const Double_t hZ = 3.;  //mm

/// error ///
Double_t sigma = 0.025;  //mm

/// config ///
TFile *ff = new TFile("boxmake.root");
TTree *g; 
Long64_t n = 0;
Int_t symnum = 0; // 0: no lock, 1: y Locked, 2: xz Locked//

Double_t val[6],val2[6],er[6],er2[6];
Double_t vstart[6],step[6],parMax[6],parMin[6],arglist[10];
Double_t Gval[6];
Double_t chi2[6]; // 0:std, 1:y axis, 2:xz axis, 3:rere, 4:min, 5: ex y

void simmig()
{ 
   TTree *a;

   TFile *fr = new TFile("results_fsimu1.root","recreate");
   TTree *tree = new TTree("tree","tree");
   Double_t dval2[6];
   tree->Branch("cx", &val2[0], "cx/D"); tree->Branch("cy", &val2[1], "cy/D"); 
   tree->Branch("cz", &val2[2], "cz/D"); tree->Branch("ta", &val2[3], "ta/D");
   tree->Branch("tb", &val2[4], "tb/D"); tree->Branch("tc", &val2[5], "tc/D");

   tree->Branch("Gx", &Gval[0], "Gx/D"); tree->Branch("Gy", &Gval[1], "Gy/D");
   tree->Branch("Gz", &Gval[2], "Gz/D"); tree->Branch("Ga", &Gval[3], "Ga/D");
   tree->Branch("Gb", &Gval[4], "Gb/D"); tree->Branch("Gc", &Gval[5], "Gc/D");

   tree->Branch("dx", &dval2[0], "dx/D"); tree->Branch("dy", &dval2[1], "dy/D");
   tree->Branch("dz", &dval2[2], "dz/D"); tree->Branch("da", &dval2[3], "da/D");
   tree->Branch("db", &dval2[4], "db/D"); tree->Branch("dc", &dval2[5], "dc/D");
   
   tree->Branch("chi2", &chi2[4], "chi2/D");

   Double_t x,y,z,x0,y0,z0;
   TVector3 v;

   TGraph2D *gg = new TGraph2D();
   Long64_t Ng = 0;

   gStyle->SetMarkerStyle(20);
   gStyle->SetMarkerSize(0.4);
   TCanvas *c1 = new TCanvas("c1","c1",1200,400);
   c1->Divide(2,1);

   for(Int_t j=0; j<1000; j++){
      ff->cd();
      g = (TTree*)ff->Get(Form("sample_%03d",j));
      n = g->GetEntries();

      a = (TTree*)ff->Get(Form("answer_%03d",j));
      a->SetBranchAddress("cx",&Gval[0]);
      a->SetBranchAddress("cy",&Gval[1]);
      a->SetBranchAddress("cz",&Gval[2]);
      a->SetBranchAddress("ta",&Gval[3]);
      a->SetBranchAddress("tb",&Gval[4]);
      a->SetBranchAddress("tc",&Gval[5]);
      a->GetEntry(0);
      
      first_fit();
      refit();

      std::cout <<"############################################" << std::endl;
      std::cout <<"############################################" << std::endl;
      std::cout << "chi^2_std(" << j << ") = " << chi2[0] << std::endl;
      std::cout <<"############################################" << std::endl;
      std::cout <<"############################################" << std::endl;
      //chi^2: 2500くらいまでいい感じ、順向き僅かに傾きで2650、7000以上はアウト、1万越えは中心も怪しい

      if(fabs(val2[4]) > 89.){
         adjfit_yaxis();
         std::cout <<"############################################" << std::endl;
         std::cout <<"############################################" << std::endl;
         std::cout << "chi^2_y(" << j << ") = " << chi2[1] << std::endl;
         std::cout <<"############################################" << std::endl;
         std::cout <<"############################################" << std::endl;
      }

      if(fabs(val2[3]) > 80. || fabs(val2[5]) > 80.){
         adjfit_xzaxis();
         std::cout <<"############################################" << std::endl;
         std::cout <<"############################################" << std::endl;
         std::cout << "chi^2_xz(" << j << ") = " << chi2[2] << std::endl;
         std::cout <<"############################################" << std::endl;
         std::cout <<"############################################" << std::endl;
      }

      if(chi2[4] > 30000.){exfit_yaxis1();}
      std::cout << "chi^2_ex1(" << j << ") = " << chi2[5] << std::endl;
      if(chi2[4] > 30000.){exfit_yaxis2();}
      std::cout << "chi^2_ex2(" << j << ") = " << chi2[5] << std::endl;

      rerefit();
      std::cout <<"############################################" << std::endl;
      std::cout <<"############################################" << std::endl;
      std::cout << "chi^2_rere(" << j << ") = " << chi2[3] << std::endl;
      std::cout << "chi^2_final(" << j << ") = " << chi2[4] << std::endl;
      std::cout <<"############################################" << std::endl;
      std::cout <<"############################################" << std::endl;

      a->GetEntry(0); //一応

      ff->cd();
      n = g->GetEntries();
      g->SetBranchAddress("x",&x0); //何故か必要
      g->SetBranchAddress("y",&y0); g->SetBranchAddress("z",&z0);
      Ng = 0;
      for(Int_t i=0; i<6; i++){
         dval2[i] = Gval[i] - val2[i];
      }
      tree->Fill();

      for(Int_t i=0; i<n; i++){
         if(i%10){continue;}
         else{
            g->GetEntry(i);
            x = x0 - val2[0]; y = y0 - val2[1]; z = z0 - val2[2];
            v = TVector3(x,y,z);
            v.RotateX(-val2[3]*DegToRad());
            v.RotateY(-val2[4]*DegToRad());
            v.RotateZ(-val2[5]*DegToRad());
            x = v.X(); y = v.Y(); z = v.Z();
            gg->SetPoint(Ng++,x,y,z);
         }
      }

      gg->SetTitle(
            Form("(%5.2lf,%5.2lf,%5.2lf)cen_(%5.2lf,%5.2lf,%5.2lf)deg; x; y; z",val2[0],val2[1],val2[2],val2[3],val2[4],val2[5])
            );
      c1->cd(1);
      g->Draw("x:y:z");
      c1->cd(2);
      gg->SetMarkerStyle(20);
      gg->SetMarkerSize(0.4);
      gg->Draw();
      c1->SaveAs(Form("./dir_simu/FinalSamples/gazou/FinalSample%03d.pdf",j));

      std::cout << "*************************************************************" << std::endl;
      std::cout << "*************************************************************" << std::endl;
      std::cout << "*************************************************************" << std::endl;
      std::cout << "------simulation " << j << " result----- " << std::endl;
      if(symnum == 0){std::cout << "Maybe No Lock" << std::endl;}
      else if(symnum == 1){std::cout << "Maybe  y Locked" << std::endl;}
      else if(symnum == 2){std::cout << "Maybe  xz Locked" << std::endl;}
      std::cout << "ParNo    :    Result    :    Answer    :    Difference" << std::endl;
      for(Int_t i=0;i<6;i++){
         std::cout << "Par_"<< i << " : " << val2[i] << " : " 
            << Gval[i] << " : " << dval2[i] <<std::endl;
      }
      std::cout << " " << std::endl;
      std::cout << "*************************************************************" << std::endl;
      std::cout << "*************************************************************" << std::endl;
      std::cout << "*************************************************************" << std::endl;
   }
   fr->cd();
   tree->Write();
   //fr->Close();
}

void first_fit()
{
   TMinuit *min = new TMinuit(6);
   min->SetPrintLevel(1);
   min->SetFCN(Loss);
   Int_t ierflg = 0;

   for(Int_t i=0;i<6;i++){
      vstart[i] = 0.; step[i] = 0.1;
      parMax[i] = 0.; parMin[i] = 0.;
   }
   for(Int_t i=3;i<6;i++){
      vstart[i] = 0.; step[i] = 0.1;
      parMax[i] = 90.; parMin[i] = -90.;
   }
   for(Int_t i=0;i<6;i++){
      min->mnparm(i, Form("val_%d",i), vstart[i], step[i], parMin[i], parMax[i], ierflg);
   }
   arglist[0] = 6.53;
   min->mnexcm("SET ERR", arglist, 1, ierflg);
   arglist[0] = 1000; 
   arglist[1] = 1;
   min->mnexcm("MINIMISE", arglist, 2, ierflg);

   for(Int_t i=0;i<6;i++){min->GetParameter(i,val[i],er[i]);}
}

void refit()
{
   TMinuit *min = new TMinuit(6);
   min->SetPrintLevel(1);
   min->SetFCN(Loss);
   Int_t ierflg = 0;

   for(Int_t i=0;i<3;i++){
      vstart[i] = val[i]; step[i] = 0.;
      parMax[i] = 0.; parMin[i] = 0.;
   }
   for(Int_t i=3;i<6;i++){
      vstart[i] = val[i]; step[i] = 0.01;
      parMax[i] = 90.; parMin[i] = -90.;
   }
   for(Int_t i=0;i<6;i++){
      min->mnparm(i, Form("val2_%d",i), vstart[i], step[i], parMin[i], parMax[i], ierflg);
   }
   arglist[0] = 3.53;
   min->mnexcm("SET ERR", arglist, 1, ierflg);
   arglist[0] = 1000; 
   arglist[1] = 1;
   min->mnexcm("MINIMISE", arglist, 2, ierflg);

   for(Int_t i=0;i<6;i++){min->GetParameter(i,val2[i],er2[i]);}
      
   Double_t fmin,fedm,errdef;
   Int_t npari,nparx,istat;
   min->mnstat(fmin,fedm,errdef,npari,nparx,istat);
   n = g->GetEntries();
   chi2[0] = fmin/(n/10-3);
   chi2[4] = chi2[0];
}

void rerefit()
{
   TMinuit *min = new TMinuit(6);
   min->SetPrintLevel(1);
   min->SetFCN(Loss);
   Int_t ierflg = 0;

   for(Int_t i=0;i<3;i++){
      vstart[i] = val2[i]; step[i] = 0.01;
      parMax[i] = 0.; parMin[i] = 0.;
   }
   for(Int_t i=3;i<6;i++){
      vstart[i] = val2[i]; step[i] = 0.;
      parMax[i] = 0.; parMin[i] = 0.;
   }
   for(Int_t i=0;i<6;i++){
      min->mnparm(i, Form("val3_%d",i), vstart[i], step[i], parMin[i], parMax[i], ierflg);
   }
   arglist[0] = 3.53;
   min->mnexcm("SET ERR", arglist, 1, ierflg);
   arglist[0] = 1000; 
   arglist[1] = 1;
   min->mnexcm("MINIMISE", arglist, 2, ierflg);

   Double_t fmin,fedm,errdef;
   Int_t npari,nparx,istat;
   min->mnstat(fmin,fedm,errdef,npari,nparx,istat);
   n = g->GetEntries();
   chi2[3] = fmin/(n/10-3);
   if(chi2[3] < chi2[4]){
      chi2[4] = chi2[3];
      for(Int_t i=0;i<6;i++){min->GetParameter(i,val2[i],er2[i]);}
   } 
}

void adjfit_yaxis()
{
   TMinuit *min = new TMinuit(6);
   min->SetPrintLevel(1);
   min->SetFCN(Loss);
   Int_t ierflg = 0;

   for(Int_t i=0;i<3;i++){
      vstart[i] = val2[i]; step[i] = 0.;
      parMax[i] = 0.; parMin[i] = 0.;
   }
   for(Int_t i=3;i<6;i++){
      vstart[i] = 0.; step[i] = 0.01;
      parMax[i] = 90.; parMin[i] = -90.;
   }
   if(val2[4] >= 80.){vstart[4] = -70.;}
   else if(val2[4] <= -80.){vstart[4] = 70.;}
      
   for(Int_t i=0;i<6;i++){
      min->mnparm(i, Form("val2y_%d",i), vstart[i], step[i], parMin[i], parMax[i], ierflg);
   }
   arglist[0] = 3.53;
   min->mnexcm("SET ERR", arglist, 1, ierflg);
   arglist[0] = 1000; 
   arglist[1] = 1;
   min->mnexcm("MINIMISE", arglist, 2, ierflg);

   Double_t fmin,fedm,errdef;
   Int_t npari,nparx,istat;
   min->mnstat(fmin,fedm,errdef,npari,nparx,istat);
   n = g->GetEntries();
   chi2[1] = fmin/(n/10-3);
   if(chi2[1] < chi2[4]){
      chi2[4] = chi2[1];
      symnum = 1;
      for(Int_t i=0;i<6;i++){min->GetParameter(i,val2[i],er2[i]);}
   }
}

void exfit_yaxis1()
{
   TMinuit *min = new TMinuit(6);
   min->SetPrintLevel(1);
   min->SetFCN(Loss);
   Int_t ierflg = 0;

   for(Int_t i=0;i<3;i++){
      vstart[i] = val2[i]; step[i] = 0.;
      parMax[i] = 0.; parMin[i] = 0.;
   }
   for(Int_t i=3;i<6;i++){
      vstart[i] = 0.;
      parMax[i] = 180.;
      parMin[i] = -180.;
   }   
   vstart[4] = 90.;
   step[3] = 0.;
   step[4] = 0.;
   step[5] = 0.01;
   for(Int_t i=0;i<6;i++){
      min->mnparm(i, Form("val2y_%d",i), vstart[i], step[i], parMin[i], parMax[i], ierflg);
   }
   arglist[0] = 1.0;
   min->mnexcm("SET ERR", arglist, 1, ierflg);
   arglist[0] = 1000; 
   arglist[1] = 1;
   min->mnexcm("MINIMISE", arglist, 2, ierflg);

   Double_t fmin,fedm,errdef;
   Int_t npari,nparx,istat;
   min->mnstat(fmin,fedm,errdef,npari,nparx,istat);
   n = g->GetEntries();

   chi2[5] = fmin/(n/10-1);
   if(chi2[5] < chi2[4]){
      chi2[4] = chi2[5];
      symnum = 1;
      for(Int_t i=0;i<6;i++){min->GetParameter(i,val2[i],er2[i]);}
   }
}

void exfit_yaxis2()
{
   TMinuit *min = new TMinuit(6);
   min->SetPrintLevel(1);
   min->SetFCN(Loss);
   Int_t ierflg = 0;

   for(Int_t i=0;i<3;i++){
      vstart[i] = val2[i]; step[i] = 0.;
      parMax[i] = 0.; parMin[i] = 0.;
   }
   for(Int_t i=3;i<6;i++){
      vstart[i] = 0.;
      parMax[i] = 180.;
      parMin[i] = -180.;
   }   
   vstart[4] = -90.;
   step[3] = 0.;
   step[4] = 0.;
   step[5] = 0.01;
   for(Int_t i=0;i<6;i++){
      min->mnparm(i, Form("val2y_%d",i), vstart[i], step[i], parMin[i], parMax[i], ierflg);
   }
   arglist[0] = 1.0;
   min->mnexcm("SET ERR", arglist, 1, ierflg);
   arglist[0] = 1000; 
   arglist[1] = 1;
   min->mnexcm("MINIMISE", arglist, 2, ierflg);

   Double_t fmin,fedm,errdef;
   Int_t npari,nparx,istat;
   min->mnstat(fmin,fedm,errdef,npari,nparx,istat);
   n = g->GetEntries();

   chi2[5] = fmin/(n/10-1);
   if(chi2[5] < chi2[4]){
      chi2[4] = chi2[5];
      symnum = 1;
      for(Int_t i=0;i<6;i++){min->GetParameter(i,val2[i],er2[i]);}
   }
}

void adjfit_xzaxis()
{
      TMinuit *min = new TMinuit(6);
      min->SetPrintLevel(1);
      min->SetFCN(Loss);
      Int_t ierflg = 0;

      for(Int_t i=0;i<3;i++){
         vstart[i] = val[i]; step[i] = 0.;
         parMax[i] = 0.; parMin[i] = 0.;
      }
      for(Int_t i=3;i<6;i++){
         vstart[i] = 0.; step[i] = 0.01;
         parMax[i] = 90.; parMin[i] = -90.;
      }

      if(val2[3] >= 80.){vstart[3] = -60.;}
      else if(val2[3] <= -80.){vstart[3] = 60.;}

      if(val2[5] >= 80.){vstart[5] = -60.;}
      else if(val2[5] <= -80.){vstart[5] = 60.;}

      for(Int_t i=0;i<6;i++){
         min->mnparm(i, Form("val2xz_%d",i), vstart[i], step[i], parMin[i], parMax[i], ierflg);
      }
      arglist[0] = 3.53;
      min->mnexcm("SET ERR", arglist, 1, ierflg);
      arglist[0] = 1000;
      arglist[1] = 1;
      min->mnexcm("MINIMISE", arglist, 2, ierflg);
      
      Double_t fmin,fedm,errdef;
      Int_t npari,nparx,istat;
      min->mnstat(fmin,fedm,errdef,npari,nparx,istat);
      n = g->GetEntries();
      chi2[2] = fmin/(n/10-3);
      if(chi2[2] < chi2[4]){
         chi2[4] = chi2[2];
         symnum = 2;
         for(Int_t i=0;i<6;i++){
            min->GetParameter(i,val2[i],er2[i]);
         }
      }
}

void Loss(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
   Double_t L = 0;
   Double_t R = 0;
   Double_t x,y,z,x0,y0,z0;
   TVector3 v;
   ff->cd();
   g->SetBranchAddress("x",&x0); //なんか必要
   g->SetBranchAddress("y",&y0); g->SetBranchAddress("z",&z0);
   n = g->GetEntries();

   for(Int_t i=0; i<n; i++){
      if(i%10){continue;}
      else{
         g->GetEntry(i);
         x = x0-par[0]; y = y0-par[1]; z = z0-par[2];
         v = TVector3(x,y,z);
         v.RotateX(-par[3]*DegToRad());
         v.RotateY(-par[4]*DegToRad());
         v.RotateZ(-par[5]*DegToRad());
         x = v.X(); y = v.Y(); z = v.Z();
         Di(x,y,z,R);
         L += (R*R)/(sigma*sigma);
      }
   }
   f = L;
}

void Di2(Double_t x, Double_t y, Double_t z, Double_t &r){
   Double_t xm,ym,zm,xp,yp,zp,zm2,zp2,zmax,zmin,zmaxr,zminr;
   Double_t R = 0.;

   xm = (x + hx) * (x + hx);
   ym = (y + hy) * (y + hy);
   zm = (z + hz) * (z + hz);
   zm2 = (z + hz2) * (z + hz2);
   xp = (x - hx) * (x - hx);
   yp = (y - hy) * (y - hy);
   zp = (z - hz) * (z - hz);
   zp2 = (z - hz2) * (z - hz2);

   zmax = -tanyz * y + (hz-hz2)/2 + 3.1;
   zmaxr = fabs(zmax - z)/sqrt(tanyz*tanyz + 1);
   zmin =  tanyz * y - (hz-hz2)/2 - 3.1;
   zminr = fabs(zmin - z)/sqrt(tanyz*tanyz + 1);

   if(y < -hy){
      if(x > hx){
         if     (z > hz) {R = sqrt(xp + ym + zp);}  //case1
         else if(z < -hz){R = sqrt(xp + ym + zm);}  //case2
         else            {R = sqrt(xp + ym);}       //case3
      }else if(x < -hx){
         if     (z > hz) {R = sqrt(xm + ym + zp);}  //case4
         else if(z < -hz){R = sqrt(xm + ym + zm);}  //case5
         else            {R = sqrt(xm + ym);}       //case6
      }else{
         if(z < -hz)     {R = sqrt(ym + zm);}       //case7
         else if(z > hz) {R = sqrt(ym + zp);}       //case8
         else            {R = sqrt(ym);}            //case9
      }
   }else if(y > hy){
      if(x > hx){
         if     (z > hz2) {R = sqrt(xp + yp + zp2);}  //case10
         else if(z < -hz2){R = sqrt(xp + yp + zm2);}  //case11
         else            {R = sqrt(xp + yp);}       //case12
      }else if(x < -hx){
         if     (z > hz2) {R = sqrt(xm + yp + zp2);}  //case13
         else if(z < -hz2){R = sqrt(xm + yp + zm2);}  //case14
         else            {R = sqrt(xm + yp);}       //case15
      }else{
         if     (z < -hz2){R = sqrt(yp + zm2);}       //case16
         else if(z > hz2) {R = sqrt(yp + zp2);}       //case17
         else            {R = sqrt(yp);}            //case18
      }
   }else{
      if(x > hx){
         if(z > zmax)    {R = sqrt(xp + zmaxr*zmaxr);}       //case19
         else if(z < zmin){R = sqrt(xp + zminr*zminr);}      //case20
         else            {R = sqrt(xp);}                     //case21
      }else if(x < -hx){
         if(z > zmax)     {R = sqrt(xm + zmaxr*zmaxr);}      //case22
         else if(z < zmin){R = sqrt(xm + zminr*zminr);}      //case23
         else            {R = sqrt(xm);}                     //case24
      }else              {R = TMath::Min(fabs(fabs(x) - hx), fabs(fabs(y)-hy));
                          R = TMath::Min(R, zminr);
                          R = TMath::Min(R, zmaxr);}         //case25
      }
   r = R;
}

void Di(Double_t x, Double_t y, Double_t z, Double_t &r){
   Double_t xm,ym,zm,xp,yp,zp;
   Double_t R[5];

   xm = (x + hX) * (x + hX);
   ym = (y + hY) * (y + hY);
   zm = (z + hZ) * (z + hZ);
   xp = (x - hX) * (x - hX);
   yp = (y - hY) * (y - hY);
   zp = (z - hZ) * (z - hZ);

   //底面
   /*
   if     (x > hX && z > hZ)              {R[5] = sqrt(xp + yp + zp);}//case1
   else if(x > hX && z < -hZ)             {R[5] = sqrt(xp + yp + zm);}//case2
   else if(x > hX && z >= -hZ && z <= hZ) {R[5] = sqrt(xp + yp);}     //case3
   else if(x < -hX && z > hZ)             {R[5] = sqrt(xm + yp + zp);}//case4
   else if(x < -hX && z < -hZ)            {R[5] = sqrt(xm + yp + zm);}//case5
   else if(x < -hX && z >= -hZ && z <= hZ){R[5] = sqrt(xm + yp);}     //case6
   else if(x >= -hX && x <= hX && z > hZ) {R[5] = sqrt(yp + zp);}     //case7
   else if(x >= -hX && x <= hX && z < -hZ){R[5] = sqrt(yp + zm);}     //case8
   else                                   {R[5] = sqrt(yp);}          //case9
   */
   //上面(USだとy<0,DSだとy>0)
   if     (x > hX && z > hZ)              {R[1] = sqrt(ym + xp + zp);}//case10
   else if(x > hX && z < -hZ)             {R[1] = sqrt(ym + xp + zm);}//case11
   else if(x > hX && z >= -hZ && z <= hZ) {R[1] = sqrt(ym + xp);}     //case12
   else if(x < -hX && z > hZ)             {R[1] = sqrt(ym + xm + zp);}//case13
   else if(x < -hX && z < -hZ)            {R[1] = sqrt(ym + xm + zm);}//case14
   else if(x < -hX && z >= -hZ && z <= hZ){R[1] = sqrt(ym + xm);}     //case15
   else if(x >= -hX && x <= hX && z < -hZ){R[1] = sqrt(ym + zm);}     //case16
   else if(x >= -hX && x <= hX && z > hZ) {R[1] = sqrt(ym + zp);}     //case17
   else                                   {R[1] = sqrt(ym);}          //case18
   //側面・小1 (x > 0)
   if     (y > hY && z > hZ)              {R[2] = sqrt(xp + yp + zp);}//case19
   else if(y > hY && z < -hZ)             {R[2] = sqrt(xp + yp + zm);}//case20
   else if(y > hY && z >= -hZ && z <= hZ) {R[2] = sqrt(xp + yp);}     //case21
   else if(y < -hY && z > hZ)             {R[2] = sqrt(xp + ym + zp);}//case22
   else if(y < -hY && z < -hZ)            {R[2] = sqrt(xp + ym + zm);}//case23
   else if(y < -hY && z >= -hZ && z <= hZ){R[2] = sqrt(xp + ym);}     //case24
   else if(y >= -hY && y <= hY && z > hZ) {R[2] = sqrt(xp + zp);}     //case25
   else if(y >= -hY && y <= hY && z < -hZ){R[2] = sqrt(xp + zm);}     //case26
   else                                   {R[2] = sqrt(xp);}          //case27
   //側面・小2 (x < 0)
   if     (y > hY && z > hZ)              {R[3] = sqrt(xm + yp + zp);}//case28
   else if(y > hY && z < -hZ)             {R[3] = sqrt(xm + yp + zm);}//case29
   else if(y > hY && z >= -hZ && z <= hZ) {R[3] = sqrt(xm + yp);}     //case30
   else if(y < -hY && z > hZ)             {R[3] = sqrt(xm + ym + zp);}//case31
   else if(y < -hY && z < -hZ)            {R[3] = sqrt(xm + ym + zm);}//case32
   else if(y < -hY && z >= -hZ && z >= hZ){R[3] = sqrt(xm + ym);}     //case33
   else if(y >= -hY && y <= hY && z > hZ) {R[3] = sqrt(xm + zp);}     //case34
   else if(y >= -hY && y <= hY && z < -hZ){R[3] = sqrt(xm + zm);}     //case35
   else                                   {R[3] = sqrt(xm);}          //case36
   //側面・大1 (z > 0)
   if     (y > hY && x > hX)              {R[4] = sqrt(xp + yp + zp);}//case37
   else if(y > hY && x < -hX)             {R[4] = sqrt(xm + yp + zp);}//case38
   else if(y > hY && x >= -hX && x <= hX) {R[4] = sqrt(yp + zp);}     //case39
   else if(y < -hY && x > hX)             {R[4] = sqrt(xp + ym + zp);}//case40
   else if(y < -hY && x < -hX)            {R[4] = sqrt(xm + ym + zp);}//case41
   else if(y < -hY && x >= -hX && x <= hX){R[4] = sqrt(ym + zp);}     //case42
   else if(y >= -hY && y <= hY && x > hX) {R[4] = sqrt(xp + zp);}     //case43
   else if(y >= -hY && y <= hY && x < -hX){R[4] = sqrt(xm + zp);}     //case44
   else                                   {R[4] = sqrt(zp);}          //case45
   //側面・大2
   if     (y > hY && x > hX)              {R[0] = sqrt(xp + yp + zm);}//case46
   else if(y > hY && x < -hX)             {R[0] = sqrt(xm + yp + zm);}//case47
   else if(y > hY && x >= -hX && x <= hX) {R[0] = sqrt(yp + zm);}     //case48
   else if(y < -hY && x > hX)             {R[0] = sqrt(xp + ym + zm);}//case49
   else if(y < -hY && x < -hX)            {R[0] = sqrt(xm + ym + zm);}//case50
   else if(y < -hY && x >= -hX && x <= hX){R[0] = sqrt(ym + zm);}     //case51
   else if(y >= -hY && y <= hY && x > hX) {R[0] = sqrt(xp + zm);}     //case52
   else if(y >= -hY && y <= hY && x < -hX){R[0] = sqrt(xm + zm);}     //case53
   else                                   {R[0] = sqrt(zm);}          //case54
   r = TMath::MinElement(5,R);
}
