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

const Double_t hx = 61.9;   //mm
const Double_t hy0 = 56.35; //mm
const Double_t hz2 = 3.1;    //mm
const Double_t hz = 5.25;  //mm
const Double_t tanyz = (Tan(ASin((hz-hz2)/hy0)));
const Double_t hy = (hz-hz2)/tanyz/2; //mm
const Double_t S1 = (2*hz2 + 2*hz) * 2*hy /2; //mm^2
const Double_t S2 = 2*hz2 * 2*hx;        //mm^2
const Double_t S3 = 2*hx * hy0;          //mm^2
const Double_t S  = 2*S1 + 2*S3 + S2;
const Double_t D2R = DegToRad();
const Double_t Pi = TMath::Pi();

Double_t x,y,z,dx,dy,dz;
TFile *f1 = new TFile("one.root");
TTree *t1 = (TTree*)f1->Get("tree1");
Long64_t N = t1->GetEntries();
TGraph2D *gr0 = new TGraph2D();
Long64_t Ng;

Double_t cx, cy, cz, cx_er, cy_er, cz_er;
Double_t tz1, ty, tz2, tz1_er, ty_er, tz2_er;
TVector3 v;

TGraph2D *g1 = new TGraph2D();

void Distance(Double_t x, Double_t y, Double_t z, Double_t &r);
void Loss(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void migrad();
void preRotate();

void init() {
   t1->SetBranchAddress("x",&x);
   t1->SetBranchAddress("y",&y);
   t1->SetBranchAddress("z",&z);
   t1->SetBranchAddress("dx",&dx);
   t1->SetBranchAddress("dy",&dy);
   t1->SetBranchAddress("dz",&dz);
}

void lack_boxmake()
{
   TFile *fout = new TFile("lack_boxmake.root","recreate");
   TTree *tr0 = new TTree("tr0","tr0");
   tr0->Branch("x",&x,"x/D");
   tr0->Branch("y",&y,"y/D");
   tr0->Branch("z",&z,"z/D");
   TTree *tr1 = new TTree("tr0","tr0");
   tr1->Branch("x",&x,"x/D");
   tr1->Branch("y",&y,"y/D");
   tr1->Branch("z",&z,"z/D");
   TCanvas *cr0 = new TCanvas("cr0", "example0",0,0,600,400);
   TGraph2D *gr0 = new TGraph2D();
   gr0->SetMarkerStyle(20);
   gr0->SetMarkerSize(0.4);
   gr0->SetNameTitle("lacking_box","lacking_box; x axis; y axis; z axis");
   Long64_t nr0 = 0;

   const Long64_t nr = 30000;
   Int_t nr1 = nr * (2*hy*2*hz)/S;
   Int_t nr2 = nr * S2/S;
   Int_t nr3 = nr * S3/S;
   
   TRandom *R = new TRandom();
   R->SetSeed();

   const Double_t E = 0.025 * 2;
   TVector3 ve;
   const Double_t p[3] = {4.48527, 0.192964, -0.00900715};

   Double_t zmax,zmin;
//// tr0 making /////
   //側面・小(x>0)
   for(Int_t i=0; i<nr1; i++){
      x = hx;
      y = -hy + 2*hy*(R->Rndm());
      if(i%10 && y < -5.){continue;}
      z = -hz + 2*hz*(R->Rndm());
      if(y > 9.){
         zmax = y*y*p[2] + y*p[1] + p[0];
      }else{zmax = 9*9*p[2] + 9*p[1] + p[0];}
      zmin = (-1) * zmax;
      if(z > zmax){continue;}
      if(z < zmin){continue;}
      else{
         v = TVector3(x,y,z);
         ve = TVector3(R->Gaus()-0.5,R->Gaus()-0.5,R->Gaus()-0.5);
         ve = E * ve;
         v = v + ve;
         x = v.X();
         y = v.Y();
         z = v.Z();
         tr0->Fill();
      }
   }
   //側面・小(x<0)
   /*for(Int_t i=0; i<nr1; i++){
      x = -hx;
      y = -hy + 2*hy*(R->Rndm());
      z = -hz + 2*hz*(R->Rndm());
      if(y > 9.){
         zmax = y*y*p[2] + y*p[1] + p[0];
      }else{zmax = 9*9*p[2] + 9*p[1] + p[0];}
      zmin = (-1) * zmax;
      if(z > zmax){continue;}                                                                       
      if(z < zmin){continue;}
      else{
         v = TVector3(x,y,z);
         ve = TVector3(R->Gaus()-0.5,R->Gaus()-0.5,R->Gaus()-0.5);
         ve = E * ve;
         v = v + ve;
         x = v.X();
         y = v.Y();
         z = v.Z();
         tr0->Fill();
      }
   }*/
   //上面
   for(Int_t i=0; i<nr2; i++){
      y = hy;
      x = -hx + 2*hx*(R->Rndm());
      z = -hz2 + 2*hz2*(R->Rndm());
      v = TVector3(x,y,z);
      ve = TVector3(R->Gaus()-0.5,R->Gaus()-0.5,R->Gaus()-0.5);
      ve = E * ve;
      v = v + ve;
      x = v.X();        
      y = v.Y();
      z = v.Z();
      tr0->Fill();
   }
   Double_t yv;
   /*//側面・大(z<0)
   for(Int_t i=0; i<nr3; i++){
      y = -hy + 2*hy*(R->Rndm());
      if(i%5 && y < -10){continue;}
      x = -hx + 2*hx*(R->Rndm());
      yv = -hy/hx * x - hy/2;
      if(y < yv){continue;}
      if(y > 9.){z = -y*y*p[2] - y*p[1] -p[0];}
      else{z = -9*9*p[2] - 9*p[1] -p[0];}
      v = TVector3(x,y,z);
      ve = TVector3(R->Gaus()-0.5,R->Gaus()-0.5,R->Gaus()-0.5);
      ve = E * ve;
      v = v + ve;
      x = v.X();
      y = v.Y();
      z = v.Z();
      tr0->Fill();
   }*/
   //側面・大(z>0)
   for(Int_t i=0; i<nr3; i++){
      y = -hy + 2*hy*(R->Rndm());
      if(i%50 && y < -5.){continue;}
      x = -hx + 2*hx*(R->Rndm());
      yv = -hy/hx * x;
      if(y > 9.){z = y*y*p[2] + y*p[1] + p[0];}
      else{z = 9*9*p[2] + 9*p[1] + p[0];}
      v = TVector3(x,y,z);
      ve = TVector3(R->Gaus()-0.5,R->Gaus()-0.5,R->Gaus()-0.5);
      ve = E * ve;
      v = v + ve;
      x = v.X();
      y = v.Y();
      z = v.Z()+(R->Rndm()-0.5)*2;
      tr0->Fill();
   }

   tr0->Write();
   Ng = tr0->GetEntries();
   for(Int_t i=0; i<Ng; i++){
      tr0->GetEntry(i);
      x *= -1;
      y *= -1;
      tr1->Fill();
      gr0->SetPoint(i,x,y,z);
   }
   gr0->Draw("P");
   gr0->Write();
   /*preRotate();
   migrad();

   Double_t Rmax = 0.;
   Double_t Rmean = 0.;
   Double_t r = 0.;

   Int_t ii =0;

   for(Int_t i=0; i<N; i++){
      if(i%100) continue;
      g1->GetPoint(i,x,y,z);
      x = x - cx;
      y = y - cy;
      z = z - cz;
      v = TVector3(x,y,z);
      v.RotateZ(tz1 * D2R);
      v.RotateY(ty * D2R);
      v.RotateZ(tz2 * D2R);
      x = v.X();
      y = v.Y();
      z = v.Z();
      g2->SetPoint(ii,x,y,z);
      ii +=1;
      Distance(x,y,z,r);
      Rmax = TMath::Max(r,Rmax);
      Rmean += r/N*100;
   }

   std::cout << "**********************************" << std::endl;
   std::cout << "Rmax = " << Rmax << std::endl;
   std::cout << "Rmean = " << Rmean << std::endl;
   std::cout << "**********************************" << std::endl;*/

///// gr0 >> gj  //////
   char* filename;
   filename = new char[200];
   TFile *f[1000];
   Double_t x0,y0,z0,C1,C2,C3,C4;
   TTree *t2[1000];
   TTree *answer[1000];
   Double_t Gx,Gy,Gz,Ga,Gb,Gc;
   //TCanvas *c2 = new TCanvas();

   for(Int_t j=0; j<1000; j++){
      Gx=(R->Rndm()-0.5)*10;
      Gy=(R->Rndm()-0.5)*10;
      Gz=(R->Rndm()-0.5)*10;
      Ga=(R->Rndm()-0.5)*10;
      Gb=(R->Rndm()-0.5)*10;
      Gc=(R->Rndm()-0.5)*10;

      //sprintf(filename, "./dir_simu/FinalSamples/FinalSample%04d.root",j);
      //f[j] = new TFile(filename,"recreate");
      t2[j] = new TTree(Form("sample_%03d",j),Form("sample_%03d",j));
      answer[j] = new TTree(Form("answer_%03d",j),Form("answer_%03d",j));
      t2[j]->Branch("x",&x);
      t2[j]->Branch("y",&y);
      t2[j]->Branch("z",&z);
      answer[j]->Branch("cx",&Gx);
      answer[j]->Branch("cy",&Gy);
      answer[j]->Branch("cz",&Gz);
      answer[j]->Branch("ta",&Ga);
      answer[j]->Branch("tb",&Gb);
      answer[j]->Branch("tc",&Gc);
      answer[j]->Fill();
      Ga = Ga * DegToRad();
      Gb = Gb * DegToRad();
      Gc = Gc * DegToRad();

      /*C1 = Cos(Ga)*Cos(Gc) - Sin(Ga)*Cos(Gb)*Sin(Gc);
      C2 = -Sin(Ga)*Cos(Gc) - Cos(Ga)*Cos(Gb)*Sin(Gc);
      C3 = Cos(Ga)*Sin(Gc) + Sin(Ga)*Cos(Gb)*Cos(Gc);
      C4 = -Sin(Ga)*Sin(Gc) + Cos(Ga)*Cos(Gb)*Cos(Gc);*/

      for(Int_t i=0; i<Ng; i++){
         tr1->GetEntry(i);
         /*x = x0*C1 + y0*C2 + z0*Sin(Gb)*Sin(Gc) + Gx;
         y = x0*C3 + y0*C4 - z0*Sin(Gb)*Cos(Gc) + Gy;
         z = x0*Sin(Ga)*Sin(Gb) + y0*Cos(Ga)*Sin(Gb) + z0*Cos(Gb) + Gz;*/
         v = TVector3(x,y,z);
         v.RotateZ(Gc);
         v.RotateY(Gb);
         v.RotateX(Ga);
         x = v.X() + Gx;
         y = v.Y() + Gy;
         z = v.Z() + Gz;
         t2[j]->Fill();
      }

      //t2[j]->Draw("x:y:z");
      fout->cd();
      t2[j]->Write();
      answer[j]->Write();
   }

   /*cr0->cd();
   gr0->Draw("P");
   cr0->SaveAs("hako.pdf");
   c2->cd();
   g2->Draw("P");
   c2->SaveAs("fitted.pdf");

   TFile *foutR0 = new TFile("boxmake.root","recreate");
   gr0->Write();
   g2->Write();
   foutR0->Close();*/
}

void Distance(Double_t x, Double_t y, Double_t z, Double_t &r){
   Double_t xi,yi,zi,d;
   Double_t D;
   for(Int_t i=0; i<Ng; i++){
      gr0->GetPoint(i,xi,yi,zi);
      d = sqrt((x-xi)*(x-xi) + (y-yi)*(y-yi) + (z-zi)*(z-zi));
      if(i==0){D = d;}
      D = TMath::Min(D, d);
   }
   r = D;
}

void Loss(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
   Double_t L = 0.;
   Double_t r = 0.;
   //TVector3 yy,zz;
   //yy = TVector3(0,1,0);
   //zz = TVector3(0,0,1);
   //yy.RotateZ(par[3] * D2R);
   //zz.Rotate(par[4] * D2R, yy);
   //sigma_i = sigma_sys*sigma_sys;

   for(Int_t i=0; i<N; i++){
      if(i%100) continue;
      g1->GetPoint(i,x,y,z);
      x = x - par[0];
      y = y - par[1];
      z = z - par[2];
      v = TVector3(x,y,z);
      v.RotateZ(par[3] * D2R);
      v.RotateY(par[4] * D2R);
      v.RotateZ(par[5] * D2R);
      x = v.X();
      y = v.Y();
      z = v.Z();
      Distance(x,y,z,r);
      //Di(x,y,z,hx,hy,hz,R);
      //sigma_i = sigma_sys*sigma_sys + dx*dx + dy*dy + dz*dz;
      //L = TMath::Max(L,r);
      L += r;
   }
   f = L;
}

void migrad()
{
   init();

   Double_t tylim = 90.;

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
   min->mnparm(4, "theta_y", vstart[4], step[4], 0, 0, ierflg);
   min->mnparm(5, "thate_z2", vstart[5], step[5], 0, 0, ierflg);

   Double_t arglist[10];
   arglist[0] = 3.53;
   min->mnexcm("SET ERR", arglist, 1, ierflg);

   arglist[0] = 1000;
   arglist[1] = 1;
   min->mnexcm("MINIMISE", arglist, 2, ierflg);

   min->GetParameter(0,cx,cx_er);
   min->GetParameter(1,cy,cy_er);
   min->GetParameter(2,cz,cz_er);
   min->GetParameter(3,tz1,tz1_er);
   min->GetParameter(4,ty,ty_er);
   min->GetParameter(5,tz2,tz2_er);
}

void preRotate()
{
   init();
   Double_t xmean = (t1->GetMaximum("x") + t1->GetMinimum("x"))/2;
   Double_t ymean = (t1->GetMaximum("y") + t1->GetMinimum("y"))/2;
   Double_t zmean = (t1->GetMaximum("z") + t1->GetMinimum("z"))/2;
   Double_t X,Y,Z;
   TVector3 vmean;
   TGraph *proj = new TGraph();
   for(Int_t i=0; i<N; i++){
      t1->GetEntry(i);
      X = x - xmean;
      Y = y - ymean;
      Z = z - zmean;
      proj->SetPoint(i,X,Y);
   }
   TF1 *line = new TF1("line", "[0]+[1]*x");
   proj->Fit(line);
   Double_t Theta[3];
   Theta[0] = TMath::ATan(line->GetParameter(1));
   for(Int_t i=0; i<N; i++){
      t1->GetEntry(i);
      v = TVector3(x-xmean,y-ymean,z-zmean);
      v.RotateZ(-Theta[0]);
      X = v.X();
      Y = v.Y();
      Z = v.Z();
      proj->SetPoint(i,Y,Z);
   }
   proj->Fit(line);
   Theta[1] = TMath::ATan(line->GetParameter(1));
   for(Int_t i=0; i<N; i++){
      t1->GetEntry(i);
      v = TVector3(x-xmean,y-ymean,z-zmean);
      v.RotateZ(-Theta[0]);
      v.RotateX(-Theta[1]);
      X = v.X();
      Y = v.Y();
      Z = v.Z();
      proj->SetPoint(i,Z,X);
   }
   proj->Fit(line);
   Theta[2] = TMath::ATan(line->GetParameter(1));
   for(Int_t i=0; i<N; i++){
      t1->GetEntry(i);
      v = TVector3(x-xmean,y-ymean,z-zmean);
      v.RotateZ(-Theta[0]);
      v.RotateX(-Theta[1]);
      v.RotateY(-Theta[2]);
      X = v.X();
      Y = v.Y();
      Z = v.Z();
      g1->SetPoint(i,X,Y,Z);
   }
}

