#include <iostream>
#include <fstream>
#include <TGraph.h>
#include <TCanvas.h>

using namespace std;

void Rndm_Graph()
{
   Double_t x,y,z;
   vector<Double_t> X,Y,Z;

   ifstream ifs("Rndm5.dat");
   while(ifs >> x >> y >> z){
      X.push_back(x);
      Y.push_back(y);
      Z.push_back(z);
   }
   ifs.close();

   Double_t r;
   TGraph *tg = new TGraph();
   for(Int_t i=0; i<512; i++){
      r = sqrt(X[i]*X[i] + Y[i]*Y[i] + Z[i]*Z[i]);
      tg->SetPoint(i,i,r);
   }

   TCanvas *c = new TCanvas("c","c",0,0,600,400);
   tg->SetMarkerStyle(8);
   tg->SetMarkerSize(0.8);
   tg->Draw("ap");
}
