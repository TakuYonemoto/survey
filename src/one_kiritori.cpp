#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
void one_kiritori()
{
   TFile *f = new TFile("one_kiritori.root");
   TTree *t1 = (TTree*)f->Get("tree");
   Double_t x, y, z;
   t1->SetBranchAddress("x",&x);
   t1->SetBranchAddress("y",&y);
   t1->SetBranchAddress("z",&z);

   TCanvas *c = new TCanvas("c", "Graph2D",0,0,600,400);
   TGraph2D *dt = new TGraph2D();
   dt->SetTitle("1ko; x axis; y axis; z axis");

   Double_t x0 = 205;
   Double_t y0 = -110;
   Long64_t N = 0;

   //read all entries and fill the histograms
   Long64_t nentries = t1->GetEntries();
   for (Long64_t i=0;i<nentries;i++) {
       t1->GetEntry(i);
       if (x > x0) continue;
       if (y < y0) continue;
       dt->SetPoint(N,x,y,z);
       N += 1;
  }
   gStyle->SetPalette(1);
   gStyle->SetMarkerStyle(20);
   dt->Draw("pcol");
   
   Double_t xmax = dt->GetXmax();
   Double_t xmin = dt->GetXmin();
   Double_t xmean = (xmax + xmin)/2;

   Double_t zmax = dt->GetZmax();
   Double_t zmin = dt->GetZmin();
   Double_t zmean = (zmax + zmin)/2;

   Double_t ymax = dt->GetYmax();
   Double_t ymin = dt->GetYmin();
   Double_t ymean = (ymax + ymin)/2;

   //std::cout << xmax << "," << ymean << "," << zmax << std::endl;

   //TF3 *plane0  = new TF3("plane0","",xmin,xmax,ymean-1,ymean+1,zmin,zmax);
   //plane0->Draw("x:y:z");

   TGraph *dt1 = new TGraph();
   dt1->SetTitle("kiritori; x axis; z axis");

   Double_t x1,y1,z1;
   Int_t N1 = 0;
   for (Int_t i=0;i<N;i++){
       dt->GetPoint(i,x1,y1,z1);
       if (y1 < ymean - 0.01) continue;
       if (y1 > ymean + 0.01) continue;
       if (z1 < 560 - 1.16 * x1) continue;
       dt1->SetPoint(N1,x1,z1);
       N1 += 1;
       //std::cout << x1 << "," << y1 << "," << z1 <<std::endl;
   }

   TCanvas *c1 = new TCanvas("c1", "Graph",0,0,600,400);
   gStyle->SetPalette(1);
   gStyle->SetMarkerStyle(20);
   dt1->Draw();
   //std::cout << "xmax=" << dt1->GetXmax() << " zmax=" << dt1->GetYmax() << std::endl;
   std::cout << "Number of Points= " << N1 << std::endl;

   TF1 *line = new TF1("line","[0]+[1]*x");
   dt1->Fit(line);
   line->Draw("same");

   Double_t p0 = line->GetParameter(0);
   Double_t p1 = line->GetParameter(1);
   Double_t p12 = p1 * p1;
   Double_t sinxz = -p1/sqrt(1+ p12);
   Double_t cosxz = -1/sqrt(1 + p12);
////////////////////////////////////////////////////////
   TGraph *dt11 = new TGraph();
   dt11->SetTitle("zenbu; x axis; z axis");

   Double_t x11,y11,z11;
   Int_t N11 = 0;
   for (Int_t i=0;i<N;i++){
       dt->GetPoint(i,x11,y11,z11);
       //Double_t xp11 = p1 * x11;
       //if (z11 < p0 + xp11 - 4.0) continue;
       dt11->SetPoint(N11,x11,z11);
       N11 += 1;
       //std::cout << x1 << "," << y1 << "," << z1 <<std::endl;                                     
   }
   TCanvas *c11 = new TCanvas("c11", "Graph",0,0,600,400);
   gStyle->SetPalette(1);
   gStyle->SetMarkerStyle(5);
   dt11->Draw();
   //std::cout << "xmax=" << dt1->GetXmax() << " zmax=" << dt1->GetYmax() << std::endl;
   std::cout << "Number of Full_fit_Points= " << N11 << std::endl;

   TF1 *line11 = new TF1("line11","[0]+[1]*x");
   line11->SetParameters(p0, p1);
   dt11->Fit(line11);
   line->Draw("same");

   Double_t p0a = line11->GetParameter(0);
   Double_t p1a = line11->GetParameter(1);
   Double_t p12a = p1a * p1a;
   Double_t sinxza = -p1a/sqrt(1+ p12a);
   Double_t cosxza = -1/sqrt(1 + p12a);
////////////////////////////////////////

   TGraph2D *dt2 = new TGraph2D();
   TGraph2D *dt22 = new TGraph2D();
   dt2->SetTitle("kaiten; x axis; y axis; z axis");
   dt22->SetTitle("kaiten_2; x axis; y axis; z axis");
   Double_t x2a,y2a,z2a,x2,y2,z2,x22,y22,z22;
   for (Int_t i=0;i<N;i++){
       dt->GetPoint(i,x2a,y2a,z2a);

       x2a = x2a - xmean;
       y2 = y2a - ymean;
       z2a = z2a - zmean;

       x2 = x2a * cosxz + z2a  * sinxz;
       z2 = x2a * (-1) * sinxz + z2a * cosxz;

       x22 = x2a * cosxza + z2a  * sinxza;
       y22 = y2a - ymean;
       z22 = x2a * (-1) * sinxza + z2a * cosxza;
       dt2->SetPoint(i,x2,y2,z2);
       dt22->SetPoint(i,x22,y22,z22);
       //std::cout << x2 << "," << y2 << "," << z2 <<std::endl;
   }
   TCanvas *c2 = new TCanvas("c2", "Graph2D",0,0,600,400);
   gStyle->SetPalette(1);
   gStyle->SetMarkerStyle(5);
   //dt2->SetMinimum(-70);
   //dt2->SetMaximum(70);
   //dt2->GetYaxis()->SetLimits(-70,70);
   dt2->Draw("pcol");

   TCanvas *c3 = new TCanvas("c3", "Graph2D",0,0,600,400);
   //dt2->SetTitle("kaiten_syukusyou; x axis; y axis; z axis");
   gStyle->SetPalette(1);
   gStyle->SetMarkerStyle(5);
   //dt22->SetMinimum(-70);
   //dt22->SetMaximum(70);
   //dt22->GetYaxis()->SetLimits(-70,70);
   dt22->Draw("pcol");

   std::cout << "Xmax=" << dt2->GetXmax() << ", Xmin=" << dt2->GetXmin() << std::endl;
   std::cout << "Ymax=" << dt2->GetYmax() << ", Ymin=" << dt2->GetYmin() << std::endl;
   std::cout << "Zmax=" << dt2->GetZmax() << ", Zmin=" << dt2->GetZmin() << std::endl;

   std::cout << "\n" << "Xmax2=" << dt22->GetXmax() << ", Xmin2=" << dt22->GetXmin() << std::endl;
   std::cout << "Ymax2=" << dt22->GetYmax() << ", Ymin2=" << dt22->GetYmin() << std::endl;
   std::cout << "Zmax2=" << dt22->GetZmax() << ", Zmin2=" << dt22->GetZmin() << std::endl;
   //TF2 *plane = new TF2("plane","[0]+[1]*x+[2]*y");
   //dt1->Fit(plane);
   //plane->Draw("surf same");

}
