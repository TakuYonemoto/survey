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
void Loss(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

const Int_t CNum = 433;
const Int_t ANum = (CNum - 256)/16;
const Int_t TNum = CNum % 2;


Double_t x,y,z;
Double_t cx[16],cy[16],cz[16],tx[16],ty[16],tz[16];
TFile *f = new TFile("prerotate.root");
TTree *t1 = (TTree*)f->Get(Form("t3_%d",CNum-256));
TTree *tree2 = (TTree*)f->Get(Form("tree%d",TNum+2));
TPolyMarker3D *gr0 = (TPolyMarker3D*)f->Get("gr0");
Long64_t N = t1->GetEntries();

Double_t D2R = TMath::DegToRad();
TFile *fo = new TFile(Form("migrad_%d.root",CNum),"recreate");
TTree *tree1 = new TTree("tree1","tree1");
TTree *tree = new TTree("par2tree","par2tree");

void init() {
t1->SetBranchAddress("x",&x);
t1->SetBranchAddress("y",&y);
t1->SetBranchAddress("z",&z);
}
void init2() {
tree2->SetBranchAddress(Form("cx%d",TNum),cx);
tree2->SetBranchAddress(Form("cy%d",TNum),cy);
tree2->SetBranchAddress(Form("cz%d",TNum),cz);
tree2->SetBranchAddress(Form("tx%d",TNum),tx);
tree2->SetBranchAddress(Form("ty%d",TNum),ty);
tree2->SetBranchAddress(Form("tz%d",TNum),tz);
tree2->GetEntry(0);
}

/// fit-box size ///
Double_t hX = 60.;
Double_t hY = 20.;
Double_t hZ = 2.5;

Double_t sigma_sys = 0.050;

void migrad_433()
{
   init();
   init2();

   TVector3 v;
   TGraph2D *ggg = new TGraph2D();
   ggg->SetNameTitle("ggg",Form("counter_%d_after_PreTrans; w[mm]; v[mm]; u[mm]",CNum));

   tree1->Branch("x",&x,"x/D");
   tree1->Branch("y",&y,"y/D");
   tree1->Branch("z",&z,"z/D");
   for(Int_t i=0; i<N; i++){
      t1->GetEntry(i);
      x = x - cx[ANum];
      y = y - cy[ANum];
      z = z - cz[ANum];
      v = TVector3(x,y,z);
      v.RotateX(tx[ANum] * D2R);
      v.RotateY(ty[ANum] * D2R);
      v.RotateZ(tz[ANum] * D2R);
      x = v.X();
      y = v.Y();
      z = v.Z();
      if(i % 3 == 0){ggg->SetPoint(i,x,y,z);}
      if(y <= -20. && fabs(ATan(z/y)) < ATan(2./19.)
           && fabs(ATan(x/y)) < ATan(60./20.)){
         tree1->Fill();
      }
      if(x <= -60. && fabs(ATan(z/x)) < ATan(2./60.)
            && fabs(ATan(y/x)) < ATan(18./60.)){
         tree1->Fill();
      }
      if(y <= 5. && fabs(ATan(x/z)) < ATan(55/2.5)
            && fabs(ATan(y/z)) < ATan(15./2.5)
            && z >= 2.3 && z <= 2.9){
         tree1->Fill();
      }
   }
   tree1->Write();

   TMinuit *minz = new TMinuit(6);
   minz->SetPrintLevel(1);
   minz->SetFCN(Loss);
   Int_t ierflg = 0;
   Double_t vstart[6];
   Double_t step[6];
   vstart[0] = 0.;
   vstart[1] = 0.;
   vstart[2] = 0.;
   vstart[3] = 0.;
   vstart[4] = 0.;
   vstart[5] = 0.;
   step[0] = 0.;
   step[1] = 0.;
   step[2] = 0.01;
   step[3] = 0.;
   step[4] = 0.;
   step[5] = 0.;
   minz->mnparm(0, "center_x", vstart[0], step[0], -10, 10, ierflg);
   minz->mnparm(1, "center_y", vstart[1], step[1], -10, 10, ierflg);
   minz->mnparm(2, "center_z", vstart[2], step[2], -10, 10, ierflg);
   minz->mnparm(3, "theta_x", vstart[3], step[3], -10, 10, ierflg);
   minz->mnparm(4, "theta_y", vstart[4], step[4], -10, 10, ierflg);
   minz->mnparm(5, "thate_z", vstart[5], step[5], -10, 10, ierflg);
   Double_t arglist[10];
   arglist[0] = 1.0;
   minz->mnexcm("SET ERR", arglist, 1, ierflg);

   arglist[0] = 1000;
   arglist[1] = 1;
   minz->mnexcm("MIGRAD", arglist, 2, ierflg);
   Double_t cx2, cy2, cz2, cx_er, cy_er, cz_er;
   Double_t tx2, ty2, tz2, tx_er, ty_er, tz_er;
   minz->GetParameter(2,cz2,cz_er);


   TMinuit *min = new TMinuit(6);
   min->SetPrintLevel(1);
   min->SetFCN(Loss);

   vstart[0] = 0.;
   vstart[1] = 0.;
   vstart[2] = cz2;
   vstart[3] = 0.;
   vstart[4] = 0.;
   vstart[5] = 0.;

   step[0] = 0.01;
   step[1] = 0.01;
   step[2] = 0.;
   step[3] = 0.01;
   step[4] = 0.01;
   step[5] = 0.01;

   min->mnparm(0, "center_x", vstart[0], step[0], -10, 10, ierflg);
   min->mnparm(1, "center_y", vstart[1], step[1], -10, 10, ierflg);
   min->mnparm(2, "center_z", vstart[2], step[2], -10, 10, ierflg);
   min->mnparm(3, "theta_x", vstart[3], step[3], -10, 10, ierflg);
   min->mnparm(4, "theta_y", vstart[4], step[4], -10, 10, ierflg);
   min->mnparm(5, "thate_z", vstart[5], step[5], -10, 10, ierflg);

   arglist[0] = 3.53;
   min->mnexcm("SET ERR", arglist, 1, ierflg);

   arglist[0] = 1000;
   arglist[1] = 1;
   min->mnexcm("MIGRAD", arglist, 2, ierflg);

   min->GetParameter(0,cx2,cx_er);
   min->GetParameter(1,cy2,cy_er);
   min->GetParameter(2,cz2,cz_er);
   min->GetParameter(3,tx2,tz_er);
   min->GetParameter(4,ty2,ty_er);
   min->GetParameter(5,tz2,tz_er);

   cout<< "center=( " << cx2 << "+/-" << cx_er << endl;
   cout<< "\t" << cy2 << "+/-" << cy_er << endl;
   cout<< "\t" << cz2 << "+/-" << cz_er << " )" <<endl;
   cout<< "thetas=(" << tx2 << "," << ty2 << "," << tz2 << ")" <<endl;

   cout<<"***********************************"<<endl;
   tree->Branch("cx",&cx[ANum],"cx/D");
   tree->Branch("cy",&cy[ANum],"cy/D");
   tree->Branch("cz",&cz[ANum],"cz/D");
   tree->Branch("tx",&tx[ANum],"tx/D");
   tree->Branch("ty",&ty[ANum],"ty/D");
   tree->Branch("tz",&tz[ANum],"tz/D");
   tree->Branch("cx2",&cx2,"cx2/D");
   tree->Branch("cy2",&cy2,"cy2/D");
   tree->Branch("cz2",&cz2,"cz2/D");
   tree->Branch("tx2",&tx2,"tx2/D");
   tree->Branch("ty2",&ty2,"ty2/D");
   tree->Branch("tz2",&tz2,"tz2/D");
   tree->Fill();
   tree->Write();

   TGraph2D *gg = new TGraph2D();
   gg->SetNameTitle("gg",Form("counter_%d_after_Fit; w[mm]; v[mm]; u[mm]",CNum));
   Long64_t Ngg = 0;

   init();
   init2();

   N = t1->GetEntries();
   for(Int_t i=0; i<N; i++){
      if(i%3) continue;
      t1->GetEntry(i);
      x = x - cx[ANum];
      y = y - cy[ANum];
      z = z - cz[ANum];
      v = TVector3(x,y,z);
      v.RotateX(tx[ANum] * D2R);
      v.RotateY(ty[ANum] * D2R);
      v.RotateZ(tz[ANum] * D2R);
      x = v.X();
      y = v.Y();
      z = v.Z();
      x = x - cx2;
      y = y - cy2;
      z = z - cz2;
      v = TVector3(x,y,z);
      v.RotateX(tx2 * D2R);
      v.RotateY(ty2 * D2R);
      v.RotateZ(tz2 * D2R);
      x = v.X();
      y = v.Y();
      z = v.Z();
      gg->SetPoint(Ngg++,x,y,z);
   }

   TPolyMarker3D *tt = new TPolyMarker3D();
   TPolyMarker3D *ttt = new TPolyMarker3D();
   tt->SetName("tt");
   ttt->SetName("ttt");
   N = tree1->GetEntries();
   for(Int_t i=0; i<N; i++){
      tree1->GetEntry(i);
      tt->SetPoint(i,x,y,z);
      x = x - cx2;
      y = y - cy2;
      z = z - cz2;
      v = TVector3(x,y,z);
      v.RotateX(tx2 * D2R);
      v.RotateY(ty2 * D2R);
      v.RotateZ(tz2 * D2R);
      x = v.X();
      y = v.Y();
      z = v.Z();
      ttt->SetPoint(i,x,y,z);
   }

   gr0->SetMarkerStyle(20);
   gr0->SetMarkerSize(0.4);
   gr0->SetMarkerColor(kMagenta-4);
   gr0->SetName("gr0");
   ggg->SetMargin(0.1);
   ggg->SetMarkerStyle(20);
   ggg->SetMarkerSize(0.4);
   ggg->SetMarkerColor(kBlack);
   tt->SetMarkerStyle(20);
   tt->SetMarkerSize(0.4);
   tt->SetMarkerColor(kGreen);
   gg->SetMarkerStyle(20);
   gg->SetMargin(0.1);
   gg->SetMarkerSize(0.4);
   gg->SetMarkerColor(kBlack);
   ttt->SetMarkerStyle(20);
   ttt->SetMarkerSize(0.4);
   ttt->SetMarkerColor(kGreen);
   gg->Write();
   tt->Write();
   ggg->Write();
   ttt->Write();
   gr0->Write();



   TCanvas *c1 = new TCanvas("c1","c1",0,0,1200,500);
   c1->Divide(2,1);
   c1->cd(1);
   ggg->Draw();
   tt->Draw("same");
   gr0->Draw("same");
   auto legend = new TLegend(0.1,0.75,0.4,0.9);
   legend->AddEntry("ggg",Form("Scanned data of counter %d",CNum),"P");
   legend->AddEntry("tt","Data points used for fitting","P");
   legend->AddEntry("gr0","Imaginary_scintilator_box","P");
   legend->Draw();
   c1->cd(2);
   gg->Draw();
   ttt->Draw("same");
   gr0->Draw("same");
   auto legend2 = new TLegend(0.1,0.75,0.4,0.9);
   legend2->AddEntry("gg",Form("Scanned data of counter %d",CNum),"P");
   legend2->AddEntry("ttt","Data points used for fitting","P");
   legend2->AddEntry("gr0","Imaginary_scintilator_box","P");
   legend2->Draw();

}

void Loss(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
   Double_t L = 0;
   Double_t R = 0;
   TVector3 v;
   Double_t sigma_i = sigma_sys*sigma_sys;
   N = tree1->GetEntries();

   for(Int_t i=0; i<N; i++){
      tree1->GetEntry(i);
      x = x - par[0];
      y = y - par[1];
      z = z - par[2];
      v = TVector3(x,y,z);
      v.RotateX(par[3] * D2R);
      v.RotateY(par[4] * D2R);
      v.RotateZ(par[5] * D2R);
      x = v.X();
      y = v.Y();
      z = v.Z();
      Di(x,y,z,R);
      //sigma_i = sigma_sys*sigma_sys + dx*dx + dy*dy + dz*dz;
      L += (R*R)/sigma_i;
   }
   f = L;
}

void Di(Double_t x, Double_t y, Double_t z, Double_t &r){                                           
   Double_t xm,ym,zm,xp,yp,zp;
   Double_t R[4];
   R[0] = 1000000000000000000.;
   R[2] = 1000000000000000000.;

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
   r = TMath::MinElement(4,R);
}
