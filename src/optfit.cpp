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
const Double_t hy = (hz-hz2)/tanyz/2; //28.154485mm
const Double_t hX = 60.;  //mm
const Double_t hY = 20.;  //mm
const Double_t hZ = 3.;  //mm

/// error ///
Double_t sigma = 0.025;  //mm

/// config ///
TFile *ff = new TFile("kiridashi.root");
TTree *g; 
Long64_t n = 0;
Int_t symnum = 0; // 0: no lock, 1: y Locked, 2: xz Locked//
Int_t n_d = 5;

Double_t val[6],val2[6],er[6],er2[6];
Double_t vstart[6],step[6],parMax[6],parMin[6],arglist[10];
Double_t chi2[6]; // 0:std, 1:y axis, 2:xz axis, 3:rere, 4:min, 5: ex y

void optfit()
{ 
   TFile *fr = new TFile("optfit.root","recreate");
   TTree *tree = new TTree("tree","tree");
   tree->Branch("cx", &val2[0], "cx/D"); tree->Branch("cy", &val2[1], "cy/D"); 
   tree->Branch("cz", &val2[2], "cz/D"); tree->Branch("ta", &val2[3], "ta/D");
   tree->Branch("tb", &val2[4], "tb/D"); tree->Branch("tc", &val2[5], "tc/D");
   tree->Branch("chi2", &chi2[4], "chi2/D");

   Double_t x,y,z,x0,y0,z0;
   TVector3 v;

   TGraph2D *gg = new TGraph2D();
   Long64_t Ng = 0;

   gStyle->SetMarkerStyle(20);
   gStyle->SetMarkerSize(0.4);
   TCanvas *c1 = new TCanvas("c1","c1",1200,400);
   c1->Divide(2,1);

   for(Int_t j=0; j<256; j++){
      ff->cd();
      g = (TTree*)ff->Get(Form("t3_%d",j));
      n = g->GetEntries();
      
      first_fit();
      refit();

      std::cout <<"############################################" << std::endl;
      std::cout <<"############################################" << std::endl;
      std::cout << "chi^2_std(" << j << ") = " << chi2[0] << std::endl;
      std::cout <<"############################################" << std::endl;
      std::cout <<"############################################" << std::endl;
      //chi^2: 25000くらいまでいい感じ、順向き僅かに傾きで26500、70000以上はアウト、10万越えは中心も怪しい

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
      tree->Fill();

      ff->cd();
      g->SetBranchAddress("x",&x0); //何故か必要
      g->SetBranchAddress("y",&y0);
      g->SetBranchAddress("z",&z0);
      Ng = 0;
      n = g->GetEntries();
      gg->Clear();

      for(Int_t i=0; i<n; i++){
         if(i%n_d){continue;}
         else{
            g->GetEntry(i);
            x = x0 - val2[0]; y = y0 - val2[1]; z = z0 - val2[2];
            v = TVector3(x,y,z);
            v.RotateX(-val2[3]*DegToRad());
            v.RotateY(-val2[4]*DegToRad());
            v.RotateZ(-val2[5]*DegToRad());
            x = v.X(); y = v.Y(); z = v.Z();
            gg->SetPoint(Ng,x,y,z);
            Ng += 1;
         }
      }

      gg->SetTitle(Form("(%5.2lf,%5.2lf,%5.2lf)cen_(%5.2lf,%5.2lf,%5.2lf)deg; x; y; z",val2[0],val2[1],val2[2],val2[3],val2[4],val2[5]));
      c1->cd(1);
      g->Draw("x:y:z");
      c1->cd(2);
      gg->SetMarkerStyle(20);
      gg->SetMarkerSize(0.4);
      gg->Draw();
      c1->SaveAs(Form("./fit/gazou%04d.pdf",j));

      std::cout << "*************************************************************" << std::endl;
      std::cout << "*************************************************************" << std::endl;
      std::cout << "*************************************************************" << std::endl;
      std::cout << "------fit " << j << " result----- " << std::endl;
      if(symnum == 0){std::cout << "Maybe No Lock" << std::endl;}
      else if(symnum == 1){std::cout << "Maybe  y Locked" << std::endl;}
      else if(symnum == 2){std::cout << "Maybe  xz Locked" << std::endl;}
      std::cout << "ParNo([0,1,2=x,y,z],[3,4,5=a,b,c])" << std::endl;
      for(Int_t i=0;i<6;i++){
         std::cout << "Par_"<< i << " = " << val2[i] <<std::endl;
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
   chi2[0] = fmin/(n/n_d-3);
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
   chi2[3] = fmin/(n/n_d-3);
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
   chi2[1] = fmin/(n/n_d-3);
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

   chi2[5] = fmin/(n/n_d-1);
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

   chi2[5] = fmin/(n/n_d-1);
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
      chi2[2] = fmin/(n/n_d-3);
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
   g->SetBranchAddress("y",&y0); 
   g->SetBranchAddress("z",&z0);
   n = g->GetEntries();

   for(Int_t i=0; i<n; i++){
      if(i%n_d){continue;}
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

const Double_t xmax_SidePlane = 18.1544845 * 2;
const Double_t xmin_SidePlane = -0.84551545; //x2
const Double_t x_TopPlane = -20.;            //x1

const Double_t y_SidePlane = 5.4923668;      //y2
const Double_t y_TopPlane = 3.1;             //y1

const Double_t z_PCB = 61.65;

const Double_t slope_approx = -0.12489852; // (y2-y1)/(x2-x1) : approximate side curve by line
const Double_t const_approx = 7.9903372;  // slope_approx * x1 + y1 : const of line
const Double_t slope_square = 0.015599641; // slope_approx^2


void Di2(Double_t x, Double_t y, Double_t z, Double_t &r){
   Double_t xm,ym,zm,xp,yp,zp,zm2,zp2,zmax,zmin,zmaxr,zminr;
   Double_t R[9];

   Double_t xS1 = Power(x - xmax_SidePlane, 2);
   Double_t xS2 = Power(x - xmin_SidePlane, 2);
   Double_t xT = Power(x - x_TopPlane, 2);

   Double_t yS1 = Power(y - y_SidePlane,2);
   Double_t yS2 = Power(y + y_SidePlane,2);
   Double_t yT = Power(Abs(y) - y_TopPlane, 2);

   Double_t zS1 = Power(z - z_PCB,2);
   Double_t zS2 = Power(z + z_PCB,2);

   Double_t x_th1 = (x + slope_approx * (y - const_approx))/(slope_square + 1);
   Double_t x_th2 = (x - slope_approx * (y + const_approx))/(slope_square + 1);
   Double_t y_th1 = slope_approx * x_th1 + const_approx;
   Double_t y_th2 = -slope_approx * x_th2 - const_approx;
   Double_t d_approx1 = Power((slope_approx * x - y + const_approx),2)/(slope_square + 1);
   Double_t d_approx2 = Power((-slope_approx * x - y - const_approx),2)/(slope_square + 1);
   Double_t y_th = slope_approx * x + const_approx;

   //平面・大1(y>0)
   if     (x >= xmax_SidePlane && Abs(z) >= z_PCB)       {R[0] = Sqrt(xS1 + yS1 + Min(zS1,zS2));}
   else if(x <= xmin_SidePlane && Abs(z) >= z_PCB)       {R[0] = Sqrt(xS2 + yS1 + Min(zS1,zS2));}
   else if(x >= xmax_SidePlane)                          {R[0] = Sqrt(xS1 + yS1);}
   else if(x <= xmin_SidePlane)                          {R[0] = Sqrt(xS2 + yS1);}
   else if(Abs(z) >= z_PCB)                              {R[0] = Sqrt(yS1 + Min(zS1,zS2));}
   else                                                  {R[0] = Sqrt(yS1);}
   //平面・大2(y<0)
   if     (x >= xmax_SidePlane && Abs(z) >= z_PCB)       {R[1] = Sqrt(xS1 + yS2 + Min(zS1,zS2));}
   else if(x <= xmin_SidePlane && Abs(z) >= z_PCB)       {R[1] = Sqrt(xS2 + yS2 + Min(zS1,zS2));}
   else if(x >= xmax_SidePlane)                          {R[1] = Sqrt(xS1 + yS2);}
   else if(x <= xmin_SidePlane)                          {R[1] = Sqrt(xS2 + yS2);}
   else if(Abs(z) >= z_PCB)                              {R[1] = Sqrt(yS2 + Min(zS1,zS2));}
   else                                                  {R[1] = Sqrt(yS2);}
   //平面・小1(z>0)
   if     (x >= xmax_SidePlane && Abs(y) >= y_SidePlane) {R[2] = Sqrt(xS1 + Min(yS1,yS2) + zS1);}
   else if(x <= xmin_SidePlane && Abs(y) >= y_SidePlane){R[2] = Sqrt(xS2 + Min(yS1,yS2) + zS1);}
   else if(x >= xmax_SidePlane)                          {R[2] = Sqrt(xS1 + zS1);}
   else if(x <= xmin_SidePlane)                          {R[2] = Sqrt(xS2 + zS1);}
   else if(Abs(y) >= y_SidePlane)                        {R[2] = Sqrt(Min(yS1,yS2) + zS1);}
   else                                                  {R[2] = Sqrt(zS1);}
   //平面・小2(z<0)
   if     (x >= xmax_SidePlane && Abs(y) >= y_SidePlane) {R[3] = Sqrt(xS1 + Min(yS1,yS2) + zS2);}
   else if(x <= xmin_SidePlane && Abs(y) >= y_SidePlane) {R[3] = Sqrt(xS2 + Min(yS1,yS2) + zS2);}
   else if(x >= xmax_SidePlane)                          {R[3] = Sqrt(xS1 + zS2);}
   else if(x <= xmin_SidePlane)                          {R[3] = Sqrt(xS2 + zS2);}
   else if(Abs(y) >= y_SidePlane)                        {R[3] = Sqrt(Min(yS1,yS2) + zS2);}
   else                                                  {R[3] = Sqrt(zS2);}
   //上面
   if     (Abs(y) > y_TopPlane && Abs(z) > z_PCB)        {R[4] = Sqrt(xT + yT + Min(zS1,zS2));}
   else if(Abs(y) > y_TopPlane)                          {R[4] = Sqrt(xT + yT);}
   else if(Abs(z) > z_PCB)                               {R[4] = Sqrt(xT + Min(zS1,zS2));}
   else                                                  {R[4] = Sqrt(xT);}
   //曲面・大1(y>0)
   if     (x_th1 >= xmin_SidePlane && Abs(z) >= z_PCB)   {R[5] = Sqrt(xS2 + yS1 + Min(zS1,zS2));}
   else if(x_th1 <= x_TopPlane && Abs(z) >= z_PCB)       {R[5] = Sqrt(xT + yT + Min(zS1,zS2));}
   else if(x_th1 >= xmin_SidePlane)                      {R[5] = Sqrt(xS2 + yS1);}
   else if(x_th1 <= x_TopPlane)                          {R[5] = Sqrt(xT + yT);}
   else if(Abs(z) >= z_PCB)                              {R[5] = Sqrt(d_approx1 + Min(zS1,zS2));}
   else                                                  {R[5] = Sqrt(d_approx1);}
   //曲面・大2(y<0)
   if     (x_th2 >= xmin_SidePlane && Abs(z) >= z_PCB)   {R[6] = Sqrt(xS2 + yS2 + Min(zS1,zS2));}
   else if(x_th2 <= x_TopPlane && Abs(z) >= z_PCB)       {R[6] = Sqrt(xT + yT + Min(zS1,zS2));}
   else if(x_th2 >= xmin_SidePlane)                      {R[6] = Sqrt(xS2 + yS2);}
   else if(x_th2 <= x_TopPlane)                          {R[6] = Sqrt(xT + yT);}
   else if(Abs(z) >= z_PCB)                              {R[6] = Sqrt(d_approx2 + Min(zS1,zS2));}
   else                                                  {R[6] = Sqrt(d_approx2);}
   //曲面・小1(z>0)
   if     (x >= xmin_SidePlane && Abs(y) >= y_SidePlane) {R[7] = Sqrt(xS2 + Min(yS1,yS2) + zS1);}
   else if(x <= x_TopPlane && Abs(y) >= y_TopPlane)      {R[7] = Sqrt(xT + yT + zS1);}
   else if(x >= xmin_SidePlane)                          {R[7] = Sqrt(xS2 + zS1);}
   else if(x <= x_TopPlane)                              {R[7] = Sqrt(xT + zS1);}
   else if(Abs(y) <= Abs(y_th))                          {R[7] = Sqrt(zS1);}
   else if(y_th1 >= y_SidePlane)                         {R[7] = Sqrt(xS2 + yS1 + zS1);}
   else if(y_th2 <= -y_SidePlane)                        {R[7] = Sqrt(xS2 + yS2 + zS1);}
   else                                                  {R[7] = Sqrt(Min(d_approx1,d_approx2) + zS1);}
   //曲面・小2(z<0)
   if     (x >= xmin_SidePlane && Abs(y) >= y_SidePlane) {R[8] = Sqrt(xS2 + Min(yS1,yS2) + zS2);}
   else if(x <= x_TopPlane && Abs(y) >= y_TopPlane)      {R[8] = Sqrt(xT + yT + zS2);}
   else if(x >= xmin_SidePlane)                          {R[8] = Sqrt(xS2 + zS2);}
   else if(x <= x_TopPlane)                              {R[8] = Sqrt(xT + zS2);}
   else if(Abs(y) <= Abs(y_th))                          {R[8] = Sqrt(zS2);}
   else if(y_th1 >= y_SidePlane)                         {R[8] = Sqrt(xS2 + yS1 + zS2);}
   else if(y_th2 <= -y_SidePlane)                        {R[8] = Sqrt(xS2 + yS2 + zS2);}
   else                                                  {R[8] = Sqrt(Min(d_approx1,d_approx2) + zS2);}
   r = MinElement(9,R);
}

void Di(Double_t x, Double_t y, Double_t z, Double_t &r){
   Double_t xm,ym,zm,xp,yp,zp;
   Double_t R[3];

   xm = (x + hX) * (x + hX);
   ym = (y + hY) * (y + hY);
   zm = (z + hZ) * (z + hZ);
   xp = (x - hX) * (x - hX);
   yp = (y - hY) * (y - hY);
   zp = (z - hZ) * (z - hZ);

   //底面(データが無い)
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
   /*if     (y > hY && z > hZ)              {R[2] = sqrt(xp + yp + zp);}//case19
   else if(y > hY && z < -hZ)             {R[2] = sqrt(xp + yp + zm);}//case20
   else if(y > hY && z >= -hZ && z <= hZ) {R[2] = sqrt(xp + yp);}     //case21
   else if(y < -hY && z > hZ)             {R[2] = sqrt(xp + ym + zp);}//case22
   else if(y < -hY && z < -hZ)            {R[2] = sqrt(xp + ym + zm);}//case23
   else if(y < -hY && z >= -hZ && z <= hZ){R[2] = sqrt(xp + ym);}     //case24
   else if(y >= -hY && y <= hY && z > hZ) {R[2] = sqrt(xp + zp);}     //case25
   else if(y >= -hY && y <= hY && z < -hZ){R[2] = sqrt(xp + zm);}     //case26
   else                                   {R[2] = sqrt(xp);}          //case27
   */
   //側面・小2 (x < 0)
   if     (y > hY && z > hZ)              {R[0] = sqrt(xm + yp + zp);}//case28
   else if(y > hY && z < -hZ)             {R[0] = sqrt(xm + yp + zm);}//case29
   else if(y > hY && z >= -hZ && z <= hZ) {R[0] = sqrt(xm + yp);}     //case30
   else if(y < -hY && z > hZ)             {R[0] = sqrt(xm + ym + zp);}//case31
   else if(y < -hY && z < -hZ)            {R[0] = sqrt(xm + ym + zm);}//case32
   else if(y < -hY && z >= -hZ && z >= hZ){R[0] = sqrt(xm + ym);}     //case33
   else if(y >= -hY && y <= hY && z > hZ) {R[0] = sqrt(xm + zp);}     //case34
   else if(y >= -hY && y <= hY && z < -hZ){R[0] = sqrt(xm + zm);}     //case35
   else                                   {R[0] = sqrt(xm);}          //case36
   //側面・大1 (z > 0)
   if     (y > hY && x > hX)              {R[2] = sqrt(xp + yp + zp);}//case37
   else if(y > hY && x < -hX)             {R[2] = sqrt(xm + yp + zp);}//case38
   else if(y > hY && x >= -hX && x <= hX) {R[2] = sqrt(yp + zp);}     //case39
   else if(y < -hY && x > hX)             {R[2] = sqrt(xp + ym + zp);}//case40
   else if(y < -hY && x < -hX)            {R[2] = sqrt(xm + ym + zp);}//case41
   else if(y < -hY && x >= -hX && x <= hX){R[2] = sqrt(ym + zp);}     //case42
   else if(y >= -hY && y <= hY && x > hX) {R[2] = sqrt(xp + zp);}     //case43
   else if(y >= -hY && y <= hY && x < -hX){R[2] = sqrt(xm + zp);}     //case44
   else                                   {R[2] = sqrt(zp);}          //case45
   //側面・大2
   /*if     (y > hY && x > hX)              {R[0] = sqrt(xp + yp + zm);}//case46
   else if(y > hY && x < -hX)             {R[0] = sqrt(xm + yp + zm);}//case47
   else if(y > hY && x >= -hX && x <= hX) {R[0] = sqrt(yp + zm);}     //case48
   else if(y < -hY && x > hX)             {R[0] = sqrt(xp + ym + zm);}//case49
   else if(y < -hY && x < -hX)            {R[0] = sqrt(xm + ym + zm);}//case50
   else if(y < -hY && x >= -hX && x <= hX){R[0] = sqrt(ym + zm);}     //case51
   else if(y >= -hY && y <= hY && x > hX) {R[0] = sqrt(xp + zm);}     //case52
   else if(y >= -hY && y <= hY && x < -hX){R[0] = sqrt(xm + zm);}     //case53
   else                                   {R[0] = sqrt(zm);}          //case54
   */
   r = TMath::MinElement(3,R);
}
