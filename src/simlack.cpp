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
void Loss2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void simfit();
void simrefit();
void simrerefit();
void preRotate();

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

/// error ///
Double_t sigma = 0.025;  //mm

/// config ///
TFile *ff = new TFile("lack_boxmake.root");
TTree *g; 
Long64_t n = 0;
Double_t Theta[3];
Double_t xmean,ymean,zmean;

Double_t cx, cy, cz, cx_er, cy_er, cz_er;
Double_t cx2, cy2, cz2, cx_er2, cy_er2, cz_er2;
Double_t ta, tb, tc, ta_er, tb_er, tc_er;
Double_t ta2, tb2, tc2, ta_er2, tb_er2, tc_er2;
Double_t Gx,Gy,Gz,Ga,Gb,Gc,C1,C2,C3,C4;

void simlack()
{ 
   TTree *a;

   TFile *fr = new TFile("results_simlack.root","recreate");
   TTree *tree = new TTree("tree","tree");
   Double_t dx, dy, dz, da, db, dc;
   Double_t dx2, dy2, dz2, da2, db2, dc2;
   tree->Branch("dx", &dx, "dx/D");
   tree->Branch("dy", &dy, "dy/D");
   tree->Branch("dz", &dz, "dz/D");
   tree->Branch("da", &da, "da/D");
   tree->Branch("db", &db, "db/D");
   tree->Branch("dc", &dc, "dc/D");
   tree->Branch("dx2", &dx2, "dx2/D");
   tree->Branch("dy2", &dy2, "dy2/D");
   tree->Branch("dz2", &dz2, "dz2/D");
   tree->Branch("da2", &da2, "da2/D");
   tree->Branch("db2", &db2, "db2/D");
   tree->Branch("dc2", &dc2, "dc2/D");


   char* sampleNo;
   sampleNo = new char[100];
   char* answerNo;
   answerNo = new char[100];
   Double_t x,y,z,x0,y0,z0;
   TVector3 v,zz,xx;

   TGraph2D *hh = new TGraph2D();
   TGraph2D *gg = new TGraph2D();

   TCanvas *c1 = new TCanvas("c1","c1",1200,400);
   c1->Divide(2,1);

   char *foutname; foutname = new char[100];
   char *hhtitle; hhtitle = new char[100];
   char *ggtitle; ggtitle = new char[100];
   for(Int_t j=0; j<1000; j++){
      ff->cd();
      sprintf(sampleNo, "sample_%03d",j);
      sprintf(answerNo, "answer_%03d",j);
      g = (TTree*)ff->Get(sampleNo);

      a = (TTree*)ff->Get(answerNo);
      a->SetBranchAddress("cx",&Gx);
      a->SetBranchAddress("cy",&Gy);
      a->SetBranchAddress("cz",&Gz);
      a->SetBranchAddress("ta",&Ga);
      a->SetBranchAddress("tb",&Gb);
      a->SetBranchAddress("tc",&Gc);

      preRotate();
      simfit();
      simrefit();

      n = g->GetEntries();
      ta2 = ta2*DegToRad();
      tb2 = tb2*DegToRad();
      tc2 = tc2*DegToRad();

      a->GetEntry(0);
      ff->cd();
      g->SetBranchAddress("x",&x0);
      g->SetBranchAddress("y",&y0);
      g->SetBranchAddress("z",&z0);
      for(Int_t i=0; i<n; i++){
         g->GetEntry(i);
         x = x0 - cx2;
         y = y0 - cy2;
         z = z0 - cz2;
         v = TVector3(x,y,z);
         v.RotateX(ta2);
         v.RotateY(tb2);
         v.RotateZ(tc2);
         x = v.X();
         y = v.Y();
         z = v.Z();
         gg->SetPoint(i,x,y,z);
      }

      ta2 = ta2*RadToDeg();
      tb2 = tb2*RadToDeg();
      tc2 = tc2*RadToDeg();

      dx = Gx - cx;
      dy = Gy - cy;
      dz = Gz - cz;
      da = Ga - ta;
      db = Gb - tb;
      dc = Gc - tc;
      dx2 = Gx - cx2;
      dy2 = Gy - cy2;
      dz2 = Gz - cz2;
      da2 = Ga - ta2;
      db2 = Gb - tb2;
      dc2 = Gc - tc2;
      sprintf(ggtitle, "(%5.2lf,%5.2lf,%5.2lf)cen_(%5.2lf,%5.2lf,%5.2lf)deg; x; y; z",cx2,cy2,cz2,ta2,tb2,tc2);
      gg->SetTitle(ggtitle);

      sprintf(foutname, "./dir_simu/LackingSamples/gazou/LackingSample%03d.pdf",j);
      c1->cd(1);
      hh->SetMarkerStyle(20);
      hh->SetMarkerSize(0.4);
      g->Draw("x:y:z");
      c1->cd(2);
      gg->SetMarkerStyle(20);
      gg->SetMarkerSize(0.4);
      gg->Draw();
      c1->SaveAs(foutname);

      fr->cd();
      dx; dy; dz; da; db; dc;
      dx2; dy2; dz2; da2; db2; dc2;
      tree->Fill();

      std::cout << "*************************************************************" << std::endl;
      std::cout << "*************************************************************" << std::endl;
      std::cout << "*************************************************************" << std::endl;
      std::cout << "simu " << j << " result_2: "<< dx2 << ", " << dy2 << ", " << dz2 << ", " << da2 << ", " << db2 << ", " << dc2 << std::endl;
      std::cout << "*************************************************************" << std::endl;
      std::cout << "*************************************************************" << std::endl;
      std::cout << "*************************************************************" << std::endl;
   }
   fr->cd();
   tree->Write();
}

void simfit()
{
      TMinuit *min = new TMinuit(6);
      min->SetPrintLevel(1);
      min->SetFCN(Loss);
      Int_t ierflg = 0;

      Double_t vstart[6];
      vstart[0] = xmean;
      vstart[1] = ymean;
      vstart[2] = zmean;
      vstart[3] = Theta[0];
      vstart[4] = Theta[1];
      vstart[5] = Theta[2];

      Double_t step[6];
      step[0] = 0.;
      step[1] = 0.;
      step[2] = 0.01;
      step[3] = 0.;
      step[4] = 0.;
      step[5] = 0.;


      min->mnparm(0, "cx", vstart[0], step[0], 0, 0, ierflg);
      min->mnparm(1, "cy", vstart[1], step[1], 0, 0, ierflg);
      min->mnparm(2, "cz", vstart[2], step[2], -10., 10., ierflg);
      min->mnparm(3, "theta_x(ta)", vstart[3], step[3], 0, 0, ierflg);
      min->mnparm(4, "theta_y(tb)", vstart[4], step[4], 0, 0, ierflg);
      min->mnparm(5, "thata_z(tc)", vstart[5], step[5], 0, 0, ierflg);

      Double_t arglist[10];
      arglist[0] = 6.5;
      min->mnexcm("SET ERR", arglist, 1, ierflg);

      arglist[0] = 1000;
      arglist[1] = 2;
      min->mnexcm("MINIMISE", arglist, 2, ierflg);

      min->GetParameter(0,cx,cx_er);
      min->GetParameter(1,cy,cy_er);
      min->GetParameter(2,cz,cz_er);
      min->GetParameter(3,ta,ta_er);
      min->GetParameter(4,tb,tb_er);
      min->GetParameter(5,tc,tc_er);
}

void simrefit()
{
      TMinuit *min = new TMinuit(6);
      min->SetPrintLevel(1);
      min->SetFCN(Loss);
      Int_t ierflg = 0;

      Double_t vstart[6];
      vstart[0] = cx;
      vstart[1] = cy;
      vstart[2] = cz;
      vstart[3] = ta;
      vstart[4] = tb;
      vstart[5] = tc;

      Double_t step[6];
      step[0] = 0.01;
      step[1] = 0.01;
      step[2] = 0.;
      step[3] = 0.01;
      step[4] = 0.01;
      step[5] = 0.01;

      Double_t tbmax = 180.;
      Double_t tbmin = 0.;

      min->mnparm(0, "cx2", vstart[0], step[0], -10., 10., ierflg);
      min->mnparm(1, "cy2", vstart[1], step[1], -10., 10., ierflg);
      min->mnparm(2, "cz2", vstart[2], step[2], -10., 10., ierflg);
      min->mnparm(3, "theta_x2(ta2)", vstart[3], step[3], -10., 10., ierflg);
      min->mnparm(4, "theta_y2(tb2)", vstart[4], step[4], -10., 10., ierflg);
      min->mnparm(5, "thata_z2(tc2)", vstart[5], step[5], -10., 10., ierflg);

      Double_t arglist[10];
      arglist[0] = 6.53;
      min->mnexcm("SET ERR", arglist, 1, ierflg);

      arglist[0] = 10000;
      arglist[1] = 2;
      min->mnexcm("MINIMISE", arglist, 2, ierflg);

      min->GetParameter(0,cx2,cx_er2);
      min->GetParameter(1,cy2,cy_er2);
      min->GetParameter(2,cz2,cz_er2);
      min->GetParameter(3,ta2,ta_er2);
      min->GetParameter(4,tb2,tb_er2);
      min->GetParameter(5,tc2,tc_er2);
}

void simrerefit()
{
      TMinuit *min = new TMinuit(6);
      min->SetPrintLevel(1);
      min->SetFCN(Loss);
      Int_t ierflg = 0;

      Double_t vstart[6];
      vstart[0] = cx2;
      vstart[1] = cy2;
      vstart[2] = cz2;
      vstart[3] = ta2;
      vstart[4] = tb2;
      vstart[5] = tc2;

      Double_t step[6];
      step[0] = 0.01;
      step[1] = 0.01;
      step[2] = 0.01;
      step[3] = 0.;
      step[4] = 0.;
      step[5] = 0.;

      Double_t tbmax = 180.;
      Double_t tbmin = 0.;

      min->mnparm(0, "cx3", vstart[0], step[0], 0, 0, ierflg);
      min->mnparm(1, "cy3", vstart[1], step[1], 0, 0, ierflg);
      min->mnparm(2, "cz3", vstart[2], step[2], 0, 0, ierflg);
      min->mnparm(3, "theta_x3(ta3)", vstart[3], step[3], 2*tbmax, tbmin, ierflg);
      min->mnparm(4, "theta_y3(tb3)", vstart[4], step[4], tbmax, tbmin, ierflg);
      min->mnparm(5, "thata_z3(tc3)", vstart[5], step[5], 2*tbmin, tbmin, ierflg);

      Double_t arglist[10];
      arglist[0] = 3.53;
      min->mnexcm("SET ERR", arglist, 1, ierflg);

      arglist[0] = 10000;
      arglist[1] = 2;
      min->mnexcm("MINIMISE", arglist, 2, ierflg);

      min->GetParameter(0,cx,cx_er);
      min->GetParameter(1,cy,cy_er);
      min->GetParameter(2,cz,cz_er);
      min->GetParameter(3,ta,ta_er);
      min->GetParameter(4,tb,tb_er);
      min->GetParameter(5,tc,tc_er);
}

void Loss(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
   Double_t L = 0;
   Double_t R = 0;
   Double_t x,y,z,x0,y0,z0;
   TVector3 v;
   ff->cd();
   g->SetBranchAddress("x",&x0);
   g->SetBranchAddress("y",&y0);
   g->SetBranchAddress("z",&z0);
   n = g->GetEntries();

   /*C1 = Cos(par[3])*Cos(par[5]) - Sin(par[3])*Cos(par[4])*Sin(par[5]);
   C2 = -Sin(par[3])*Cos(par[5]) - Cos(par[3])*Cos(par[4])*Sin(par[5]);
   C3 = Cos(par[3])*Sin(par[5]) + Sin(par[3])*Cos(par[4])*Cos(par[5]);
   C4 = -Sin(par[3])*Sin(par[5]) + Cos(par[3])*Cos(par[4])*Cos(par[5]);*/

   for(Int_t i=0; i<n; i++){
      g->GetEntry(i);
      x = x0 - par[0];
      y = y0 - par[1];
      z = z0 - par[2];
      //x = x0*C1 + y0*C3 + z0*Sin(par[4])*Sin(par[3]);
      //y = x0*C2 + y0*C4 + z0*Sin(par[4])*Cos(par[3]);
      //z = x0*Sin(par[5])*Sin(par[4]) - y0*Cos(par[5])*Sin(par[4]) + z0*Cos(par[4]);
      v = TVector3(x,y,z);
      v.RotateX(-par[3]*DegToRad());
      v.RotateY(-par[4]*DegToRad());
      v.RotateZ(-par[5]*DegToRad());
      x = v.X();
      y = v.Y();
      z = v.Z();
      Di(x,y,z,R);
      L += (R*R)/(sigma*sigma);
   }
   f = L;
}

void Loss2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
   Double_t L = 0;
   Double_t R = 0;
   Double_t x,y,z;
   TVector3 v;
   n = g->GetEntries();

   for(int i=0; i<n; i++){
      g->GetEntry(i);
      x = x - par[0];
      y = y - par[1];
      z = z - par[2];
      v = TVector3(x,y,z);
      v.RotateZ(par[3] * TMath::DegToRad());
      v.RotateX(par[4] * TMath::DegToRad());
      v.RotateZ(par[5] * TMath::DegToRad());
      x = v.X();
      y = v.Y();
      z = v.Z();
      Di2(x,y,z,R);
      L += (R*R)/(sigma*sigma);
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
   Double_t R[6];
   R[0] = 1000000000000;
   R[3] = 1000000000000;
   R[5] = 1000000000000;

   xm = (x + hX) * (x + hX);
   ym = (y + hY) * (y + hY);
   zm = (z + hZ) * (z + hZ);
   xp = (x - hX) * (x - hX);
   yp = (y - hY) * (y - hY);
   zp = (z - hZ) * (z - hZ);

   /*//底面
   if(x > hX){
         if     (z > hZ) {R[0] = sqrt(xp + ym + zp);}  //case1
         else if(z < -hZ){R[0] = sqrt(xp + ym + zm);}  //case2
         else            {R[0] = sqrt(xp + ym);}       //case3
   }else if(x < -hX){
         if     (z > hZ) {R[0] = sqrt(xm + ym + zp);}  //case4
         else if(z < -hZ){R[0] = sqrt(xm + ym + zm);}  //case5
         else            {R[0] = sqrt(xm + ym);}       //case6
   }else{ 
         if     (z > hZ) {R[0] = sqrt(ym + zp);}       //case7
         else if(z < -hZ){R[0] = sqrt(ym + zm);}       //case8
         //if(z < -hZ)     {R = sqrt(ym + zm);}
         //else if(z > hZ) {R = sqrt(ym + zp);} 
         //else            {R = sqrt(ym);}            
   }*/
   //上面
   if(x > hX){
         if     (z > hZ) {R[1] = sqrt(xp + yp + zp);}  //case9
         else if(z < -hZ){R[1] = sqrt(xp + yp + zm);}  //case10
         else            {R[1] = sqrt(xp + yp);}       //case11
   }else if(x < -hX){
         if     (z > hZ) {R[1] = sqrt(xm + yp + zp);}  //case12
         else if(z < -hZ){R[1] = sqrt(xm + yp + zm);}  //case13
         else            {R[1] = sqrt(xm + yp);}       //case14
   }else{
         if     (z < -hZ){R[1] = sqrt(yp + zm);}       //case15
         else if(z > hZ) {R[1] = sqrt(yp + zp);}       //case16
         else            {R[1] = sqrt(yp);}            //case17
   }
   //側面・小1
   if(y > hY){
         if(z > hZ)      {R[2] = sqrt(xp + yp + zp);}       //case18
         else if(z < -hZ){R[2] = sqrt(xp + yp + zm);}       //case19
         else            {R[2] = sqrt(xp + yp);}            //case20
   }else if(y < -hY){
         if(z > hZ)      {R[2] = sqrt(xp + ym + zp);}       //case21
         else if(z < -hZ){R[2] = sqrt(xp + ym + zm);}       //case22
         else            {R[2] = sqrt(xp + ym);}            //case23
   }else{
         if(z > hZ)      {R[2] = sqrt(xp + zp);}       //case24
         else if(z < -hZ){R[2] = sqrt(xp + zm);}       //case25
         else            {R[2] = sqrt(xp);}            //case26
   }
   /*//側面・小2
   if(y > hY){
         if(z > hZ)      {R[3] = sqrt(xm + yp + zp);}       //case27
         else if(z < -hZ){R[3] = sqrt(xm + yp + zm);}       //case28
         else            {R[3] = sqrt(xm + yp);}            //case29
   }else if(y < -hY){
         if(z > hZ)      {R[3] = sqrt(xm + ym + zp);}       //case30
         else if(z < -hZ){R[3] = sqrt(xm + ym + zm);}       //case31
         else            {R[3] = sqrt(xm + ym);}            //case32
   }else{
         if(z > hZ)      {R[3] = sqrt(xm + zp);}       //case33
         else if(z < -hZ){R[3] = sqrt(xm + zm);}       //case34
         else            {R[3] = sqrt(xm);}            //case35
   }*/
   //側面・大1
   if(y > hY){
         if(x > hX)      {R[4] = sqrt(xp + yp + zp);}       //case36
         else if(x < -hX){R[4] = sqrt(xm + yp + zp);}       //case37
         else            {R[4] = sqrt(yp + zp);}            //case38
   }else if(y < -hY){
         if(x > hX)      {R[4] = sqrt(xp + ym + zp);}       //case39
         else if(x < -hX){R[4] = sqrt(xm + ym + zp);}       //case40
         else            {R[4] = sqrt(ym + zp);}            //case41
   }else{
         if(x > hX)      {R[4] = sqrt(xp + zp);}       //case42
         else if(x < -hX){R[4] = sqrt(xm + zp);}       //case43
         else            {R[4] = sqrt(zp);}            //case44
   }
   /*//側面・大2
   if(y > hY){
         if(x > hX)      {R[5] = sqrt(xp + yp + zm);}       //case45
         else if(x < -hX){R[5] = sqrt(xm + yp + zm);}       //case46
         else            {R[5] = sqrt(yp + zm);}            //case47
   }else if(y < -hY){
         if(x > hX)      {R[5] = sqrt(xp + ym + zm);}       //case48
         else if(x < -hX){R[5] = sqrt(xm + ym + zm);}       //case49
         else            {R[5] = sqrt(ym + zm);}            //case50
   }else{
         if(x > hX)      {R[5] = sqrt(xp + zm);}       //case51
         else if(x < -hX){R[5] = sqrt(xm + zm);}       //case52
         else            {R[5] = sqrt(zm);}            //case53
   }*/
   r = TMath::MinElement(6,R);
}

void preRotate()
{
   Double_t x,y,z,x0,y0,z0;
   TVector3 v;
   ff->cd();
   g->SetBranchAddress("x",&x);
   g->SetBranchAddress("y",&y);
   g->SetBranchAddress("z",&z);
   Long64_t N = g->GetEntries();

   xmean = (g->GetMaximum("x") + g->GetMinimum("x"))/2;
   ymean = (g->GetMaximum("y") + g->GetMinimum("y"))/2;
   zmean = (g->GetMaximum("z") + g->GetMinimum("z"))/2;
   Double_t X,Y,Z;
   TGraph *proj = new TGraph();
   for(Int_t i=0; i<N; i++){
      g->GetEntry(i);
      X = x - xmean;
      Y = y - ymean;
      Z = z - zmean;
      proj->SetPoint(i,Y,Z);
   }
   TF1 *line = new TF1("line", "[0]+[1]*x");
   proj->Fit(line);
   Theta[0] = TMath::ATan(line->GetParameter(1));
   for(Int_t i=0; i<N; i++){
      g->GetEntry(i);
      v = TVector3(x-xmean,y-ymean,z-zmean);
      v.RotateZ(-Theta[0]);
      X = v.X();
      Y = v.Y();
      Z = v.Z();
      proj->SetPoint(i,Z,X);
   }
   proj->Fit(line);
   Theta[1] = TMath::ATan(line->GetParameter(1));
   for(Int_t i=0; i<N; i++){
      g->GetEntry(i);
      v = TVector3(x-xmean,y-ymean,z-zmean);
      v.RotateZ(-Theta[0]);
      v.RotateX(-Theta[1]);
      X = v.X();
      Y = v.Y();
      Z = v.Z();
      proj->SetPoint(i,X,Y);
   }
   proj->Fit(line);
   Theta[2] = TMath::ATan(line->GetParameter(1));
   /*for(Int_t i=0; i<N; i++){
      g->GetEntry(i);
      v = TVector3(x-xmean,y-ymean,z-zmean);
      v.RotateZ(-Theta[0]);
      v.RotateX(-Theta[1]);
      v.RotateY(-Theta[2]);
      X = v.X();
      Y = v.Y();
      Z = v.Z();
      g1->SetPoint(i,X,Y,Z);
   }*/
}

