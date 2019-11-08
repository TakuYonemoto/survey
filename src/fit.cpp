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
const Double_t hZ = 3.;  //mm

/// error ///
Double_t sigma = 0.025;  //mm

/// config ///
TTree *g[1000];
TFile *fk;
Int_t j=0;
TTree *t[256];
Double_t xi,yi,zi;
Long64_t n = 0;
char *treename;

TMatrixD Ra(3,3),Rb(3,3),Rc(3,3),Rot(3,3);
TArrayD ra(9),rb(9),rc(9);
Double_t cx, cy, cz, cx_er, cy_er, cz_er;
Double_t cx2, cy2, cz2, cx_er2, cy_er2, cz_er2;
Double_t ta, tb, tc, ta_er, tb_er, tc_er;
Double_t ta2, tb2, tc2, ta_er2, tb_er2, tc_er2;
Double_t Gx,Gy,Gz,Ga,Gb,Gc,C1,C2,C3,C4;

void Di(Double_t x, Double_t y, Double_t z, Double_t &r);
void Di2(Double_t x, Double_t y, Double_t z, Double_t &r);
void Loss(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void Loss2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void simfit();
void simrefit();

void fit()
{ 
   Double_t x,y,z,x0,y0,z0;
   TFile *fr = new TFile("/Users/taku/Downloads/survey_data/fit.root","recreate");
   TTree *tree = new TTree("tree","tree");
   Double_t dx, dy, dz;
   tree->Branch("x", &dx, "x/D");
   tree->Branch("y", &dy, "y/D");
   tree->Branch("z", &dz, "z/D");
   Double_t da, db, dc;
   tree->Branch("a", &da, "a/D");
   tree->Branch("b", &db, "b/D");
   tree->Branch("c", &dc, "c/D");

   TTree *h[256];
   for(Int_t i=0; i<256; i++){
      h[i] = new TTree(Form("h_%d",i),Form("h_%d",i));
      h[i]->Branch("x", &x, "x/D");
      h[i]->Branch("y", &y, "y/D");
      h[i]->Branch("z", &z, "z/D");
      h[i]->Branch("x0", &x0, "x0/D");
      h[i]->Branch("y0", &y0, "y0/D");
      h[i]->Branch("z0", &z0, "z0/D");
   }

   char *foutname;
   foutname = new char[100];
   TVector3 v,cv,Gv,dv,zz,xx;

   TCanvas *c1 = new TCanvas("c1","c1",1200,400);
   c1->Divide(2,1);

   treename = new char[100];
   fk = new TFile("kiridashi.root");

   for(j=0; j<256; j++){
      sprintf(foutname, "./fit/gazou%04d.pdf",j);

      sprintf(treename, "t3_%d", j);
      t[j] = (TTree*)fk->Get(treename);
      simfit();
      std::cout << "n_" << j << " = " << n << std::endl;

      t[j]->SetBranchAddress("x",&x0);
      t[j]->SetBranchAddress("y",&y0);
      t[j]->SetBranchAddress("z",&z0);
      zz = TVector3(0,0,1);
      xx = TVector3(1,0,0);
      xx.RotateZ(ta);
      zz.Rotate(tb,xx);
      /*ra.Reset(0);
      rb.Reset(0);
      rc.Reset(0);
      ra[0] = Cos(ta);
      ra[1] = -Sin(ta);
      ra[3] = Sin(ta);
      ra[4] = Cos(ta);
      ra[8] = 1.;
      Ra.SetMatrixArray(ra.GetArray());
      rb[4] = Cos(tb);
      rb[8] = Cos(tb);
      rb[5] = -Sin(tb);
      rb[7] = Sin(tb);
      rb[0] = 1.;
      Rb.SetMatrixArray(rb.GetArray());
      rc[0] = Cos(tc);
      rc[1] = -Sin(tc);
      rc[3] = Sin(tc);
      rc[4] = Cos(tc);
      rc[8] = 1.;
      Rc.SetMatrixArray(rc.GetArray());
      Rot = Rb*Ra;
      Rot = Rc*Rot;*/
      /*C1 = Cos(ta)*Cos(tc) - Sin(ta)*Cos(tb)*Sin(tc);
      C2 = -Sin(ta)*Cos(tc) - Cos(ta)*Cos(tb)*Sin(tc);
      C3 = Cos(ta)*Sin(tc) + Sin(ta)*Cos(tb)*Cos(tc);
      C4 = -Sin(ta)*Sin(tc) + Cos(ta)*Cos(tb)*Cos(tc);*/

      for(Int_t i=0; i<n; i++){
         t[j]->GetEntry(i);
         x0 = x0 - cx;
         y0 = y0 - cy;
         z0 = z0 - cz;
         v = TVector3(x0,y0,z0);
         v.RotateZ(ta);
         v.Rotate(tb,xx);
         v.Rotate(tc,zz);
         //v = Rot * v;
         /*v.RotateZ(ta);
         v.RotateX(tb);
         v.RotateZ(tc);*/
         x = v.X();
         y = v.Y();
         z = v.Z();
         /*x = C1*x0 + C3*y0 + z0*Sin(tb)*Sin(ta);
         y = C2*x0 + C4*y0 + z0*Sin(tb)*Cos(ta);
         z = Sin(tc)*Sin(tb)*x0 - Cos(tc)*Sin(tb)*y0 + Cos(tb)*z0;*/
         h[j]->Fill();
      }

      c1->cd(1);
      t[j]->Draw("x:y:z");
      c1->cd(2);
      h[j]->Draw("x:y:z");
      c1->SaveAs(foutname);

      ta = ta*RadToDeg();
      tb = tb*RadToDeg();
      tc = tc*RadToDeg();

      dx = cx;
      dy = cy;
      dz = cz;
      da = ta;
      db = tb;
      dc = tc;

      fr->cd();
      dx; dy; dz; da; db; dc;
      tree->Fill();

      std::cout << "*************************************************************" << std::endl;
      std::cout << "*************************************************************" << std::endl;
      std::cout << "*************************************************************" << std::endl;
      std::cout << "result " << j << ": " << dx << ", " << dy << ", " << dz << ", " << da << ", " << db << ", " << dc << std::endl;
      std::cout << "*************************************************************" << std::endl;
      std::cout << "*************************************************************" << std::endl;
      std::cout << "*************************************************************" << std::endl;
   }
   fr->cd();
   tree->Write();
   for(Int_t i=0; i<256; i++){
      h[i]->Write();
   }
   fr->Close();
}

void simfit()
{
      TMinuit *min = new TMinuit(6);
      min->SetPrintLevel(1);
      min->SetFCN(Loss);
      Int_t ierflg = 0;

      Double_t vstart[6];
      vstart[0] = 0.;
      vstart[1] = 0.;
      vstart[2] = 0.;
      vstart[3] = 0.;
      vstart[4] = 0.;
      vstart[5] = 0.;

      Double_t step[6];
      step[0] = 0.1;
      step[1] = 0.1;
      step[2] = 0.1;
      step[3] = 0.1;
      step[4] = 0.1;
      step[5] = 0.1;

      min->mnparm(0, "x0", vstart[0], step[0], 0, 0, ierflg);
      min->mnparm(1, "y0", vstart[1], step[1], 0, 0, ierflg);
      min->mnparm(2, "z0", vstart[2], step[2], 0, 0, ierflg);
      min->mnparm(3, "theta_z1", vstart[3], step[3], -TMath::Pi(), TMath::Pi(), ierflg);
      min->mnparm(4, "theta_x", vstart[4], step[4], -TMath::Pi(), TMath::Pi(), ierflg);
      min->mnparm(5, "thate_z2", vstart[5], step[5], -TMath::Pi(), TMath::Pi(), ierflg);

      Double_t arglist[10];
      arglist[0] = 3.53;
      min->mnexcm("SET ERR", arglist, 1, ierflg);

      arglist[0] = 1000;
      arglist[1] = 1;
      min->mnexcm("MIGRAD", arglist, 2, ierflg);

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
      min->SetFCN(Loss2);
      Int_t ierflg = 0;

      Double_t vstart[6];
      vstart[0] = 0.;
      vstart[1] = 0.;
      vstart[2] = 0.;
      vstart[3] = 0.;
      vstart[4] = 0.;
      vstart[5] = 0.;

      Double_t step[6];
      step[0] = 0.01;
      step[1] = 0.01;
      step[2] = 0.01;
      step[3] = 0.01;
      step[4] = 0.01;
      step[5] = 0.01;

      min->mnparm(0, "x0", vstart[0], step[0], 0, 0, ierflg);
      min->mnparm(1, "y0", vstart[1], step[1], 0, 0, ierflg);
      min->mnparm(2, "z0", vstart[2], step[2], 0, 0, ierflg);
      min->mnparm(3, "theta_z1", vstart[3], step[3], 0, 0, ierflg);
      min->mnparm(4, "theta_x", vstart[4], step[4], 0, 0, ierflg);
      min->mnparm(5, "thate_z2", vstart[5], step[5], 0, 0, ierflg);

      Double_t arglist[10];
      arglist[0] = 3.53;
      min->mnexcm("SET ERR", arglist, 1, ierflg);

      arglist[0] = 1000;
      arglist[1] = 1;
      min->mnexcm("MINIMISE", arglist, 2, ierflg);

      min->GetParameter(0,cx2,cx_er2);
      min->GetParameter(1,cy2,cy_er2);
      min->GetParameter(2,cz2,cz_er2);
      min->GetParameter(3,ta2,ta_er2);
      min->GetParameter(4,tb2,tb_er2);
      min->GetParameter(5,tc2,tc_er2);
}


void Loss(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
   Double_t L = 0;
   Double_t R = 0;
   Double_t x1,y1,z1,x10,y10,z10;
   TVector3 v,xx,zz;
   t[j]->SetBranchAddress("x",&x10);
   t[j]->SetBranchAddress("y",&y10);
   t[j]->SetBranchAddress("z",&z10);
   n = t[j]->GetEntries();
   xx = TVector3(1,0,0);
   zz = TVector3(0,0,1);
   xx.RotateZ(par[3]);
   zz.Rotate(par[4],xx);
   /*ra.Reset(0);
   rb.Reset(0);
   rc.Reset(0);
   ra[0] = Cos(par[3]);
   ra[1] = -Sin(par[3]);
   ra[3] = Sin(par[3]);
   ra[4] = Cos(par[3]);
   ra[8] = 1.;
   Ra.SetMatrixArray(ra.GetArray());
   rb[4] = Cos(par[4]);
   rb[8] = Cos(par[4]);
   rb[5] = -Sin(par[4]);
   rb[7] = Sin(par[4]);
   rb[0] = 1.;
   Rb.SetMatrixArray(rb.GetArray());
   rc[0] = Cos(par[5]);
   rc[1] = -Sin(par[5]);
   rc[3] = Sin(par[5]);
   rc[4] = Cos(par[5]);
   rc[8] = 1.;
   Rc.SetMatrixArray(rc.GetArray());
   Rot = Rb*Ra;
   Rot = Rc*Rot;*/
   /*C1 = Cos(par[3])*Cos(par[5]) - Sin(par[3])*Cos(par[4])*Sin(par[5]);
   C2 = -Sin(par[3])*Cos(par[5]) - Cos(par[3])*Cos(par[4])*Sin(par[5]);
   C3 = Cos(par[3])*Sin(par[5]) + Sin(par[3])*Cos(par[4])*Cos(par[5]);
   C4 = -Sin(par[3])*Sin(par[5]) + Cos(par[3])*Cos(par[4])*Cos(par[5]);*/

   for(Int_t i=0; i<n; i++){
      t[j]->GetEntry(i);
      x10 = x10 - par[0];
      y10 = y10 - par[1];
      z10 = z10 - par[2];
      v = TVector3(x10,y10,z10);
      /*v.RotateZ(par[3]);
      v.RotateX(par[4]);
      v.RotateZ(par[5]);*/
      //v = Rot * v;
      v.RotateZ(par[3]);
      v.Rotate(par[4],xx);
      v.Rotate(par[5],zz);
      x1 = v.X();
      y1 = v.Y();
      z1 = v.Z();
      /*x = x0*C1 + y0*C3 + z0*Sin(par[4])*Sin(par[3]);
      y = x0*C2 + y0*C4 + z0*Sin(par[4])*Cos(par[3]);
      z = x0*Sin(par[5])*Sin(par[4]) - y0*Cos(par[5])*Sin(par[4]) + z0*Cos(par[4]);*/
      Di(x1,y1,z1,R);
      L += (R*R)/(sigma*sigma);
   }
   f = L;
}

void Loss2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
   Double_t L = 0;
   Double_t R = 0;
   Double_t x,y,z;
   TVector3 v;
   n = t[j]->GetEntries();

   for(Int_t i=0; i<n; i++){
      t[j]->GetEntry(i);
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
   Double_t R = 0.;

   xm = (x + hX) * (x + hX);
   ym = (y + hY) * (y + hY);
   zm = (z + hZ) * (z + hZ);
   xp = (x - hX) * (x - hX);
   yp = (y - hY) * (y - hY);
   zp = (z - hZ) * (z - hZ);

   if(y < -hY){
      if(x > hX){
         if     (z > hZ) {R = sqrt(xp + ym + zp);}  //case1
         else if(z < -hZ){R = sqrt(xp + ym + zm);}  //case2
         else            {R = sqrt(xp + ym);}       //case3
      }else if(x < -hX){
         if     (z > hZ) {R = sqrt(xm + ym + zp);}  //case4
         else if(z < -hZ){R = sqrt(xm + ym + zm);}  //case5
         else            {R = sqrt(xm + ym);}       //case6
      }else{
         if(z < -hZ)     {R = sqrt(ym + zm);}       //case7
         else if(z > hZ) {R = sqrt(ym + zp);}       //case8
         else            {R = sqrt(ym);}            //case9
      }
   }else if(y > hY){
      if(x > hX){
         if     (z > hZ) {R = sqrt(xp + yp + zp);}  //case10
         else if(z < -hZ){R = sqrt(xp + yp + zm);}  //case11
         else            {R = sqrt(xp + yp);}       //case12
      }else if(x < -hX){
         if     (z > hZ) {R = sqrt(xm + yp + zp);}  //case13
         else if(z < -hZ){R = sqrt(xm + yp + zm);}  //case14
         else            {R = sqrt(xm + yp);}       //case15
      }else{
         if     (z < -hZ){R = sqrt(yp + zm);}       //case16
         else if(z > hZ) {R = sqrt(yp + zp);}       //case17
         else            {R = sqrt(yp);}            //case18
      }
   }else{
      if(x > hX){
         if(z > hZ)      {R = sqrt(xp + zp);}       //case19
         else if(z < -hZ){R = sqrt(xp + zm);}       //case20
         else            {R = sqrt(xp);}            //case21
      }else if(x < -hX){
         if(z > hZ)      {R = sqrt(xm + zp);}       //case22
         else if(z < -hZ){R = sqrt(xm + zm);}       //case23
         else            {R = sqrt(xm);}            //case24
      }else              {R = TMath::Min(fabs(fabs(z) - hZ),fabs(fabs(x) - hX)); //case25
                          R = TMath::Min(fabs(fabs(y) - hY), R);
      }
   }
   r = R;
}
