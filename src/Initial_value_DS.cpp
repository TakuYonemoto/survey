#include <TMath.h>
#include <TVector3.h>
#include <TTree.h>
#include <TFile.h>
#include <TPolyMarker3D.h>
#include <TGraph2D.h>
#include <TRandom.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>

using namespace std;

void Initial_value_DS()
{
   vector<TVector3> v_geo;
   vector<Double_t> phi_geo, theta_geo, psi_geo;
   vector<Int_t> id_geo,idx_geo,order_geo,hitorder_geo;
   const Int_t bufsize = 256;

   Int_t id,idx,order,hitorder;
   Double_t x,y,z,phi,theta,psi;
   TVector3 v;
   string str;
   ifstream ifs("data/Geometry.csv");
   if(!ifs){
      cout << "Error: File cannot open!" << endl;
   } 
   while(getline(ifs,str)){
      sscanf(str.data(), "%d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%d,%d",&id,&idx,&x,&y,&z,&phi,&theta,&psi,&order,&hitorder);
      x *= 10.;
      y *= 10.;
      z *= 10.;
      v = TVector3(x,y,z);
      v_geo.push_back(v);
      phi_geo.push_back(phi);
      theta_geo.push_back(theta);
      psi_geo.push_back(psi);
      id_geo.push_back(id);
      idx_geo.push_back(idx);
      order_geo.push_back(order);
      hitorder_geo.push_back(hitorder);
   }
   ifs.close();

   vector<Double_t> par, par_err;
   ifs.open("outputfiles/FARO_COBRA_calib2018_DS.csv");
   if(!ifs){
      cout << "Error: File cannot open!" << endl;
   } 
   while(getline(ifs,str)){
      sscanf(str.data(),"%lf,%lf",&x,&y);
      par.push_back(x);
      par_err.push_back(y);
   }
   ifs.close();
   const TVector3 v_par(par[0],par[1],par[2]);

   TRandom *generator = new TRandom();
   generator->SetSeed();
   
   TPolyMarker3D *pm3 = new TPolyMarker3D();
   pm3->SetName("pm3origin");
   pm3->SetMarkerStyle(8);
   pm3->SetMarkerSize(0.15);
   pm3->SetMarkerColor(kBlue);
   TPolyMarker3D *pm3h = new TPolyMarker3D();
   pm3h->SetName("pm3origin_5cm");
   pm3h->SetMarkerStyle(8);
   pm3h->SetMarkerSize(0.15);
   pm3h->SetMarkerColor(kBlue);
   TFile *boxmake1 = new TFile("outputfiles/boxmake1.root","read");
   TPolyMarker3D *ppm3 = (TPolyMarker3D*)boxmake1->Get("pm3");
   Int_t N = ppm3->GetN();
   for(Int_t i=0; i < N; i++){
      ppm3->GetPoint(i,x,y,z);
      pm3->SetPoint(i,y+8.154485,z,x);
      pm3h->SetPoint(i,y+8.154485-10.,z,x);
   }
   for(Int_t i=0; i < 3000; i++){
      x = 0;
      y = -3.1 + (generator->Rndm()*6.2);
      z = -60. + (generator->Rndm()*120.);
      pm3->SetPoint(N++,x,y,z);
   }
   for(Int_t i=0; i < 3000; i++){
      x = -20. + (generator->Rndm()*40.);
      y = -3.1 + (generator->Rndm()*6.2);
      z = 0;
      pm3->SetPoint(N++,x,y,z);
   }
   for(Int_t i=0; i < 3000; i++){
      x = -20. + (generator->Rndm()*40.);
      y = 0;
      z = -60. + (generator->Rndm()*120.);
      pm3->SetPoint(N++,x,y,z);
   }
   N = ppm3->GetN();
   for(Int_t i=0; i < 3000; i++){
      x = 0;
      y = -3.1 + (generator->Rndm()*6.2);
      z = -60. + (generator->Rndm()*120.);
      pm3h->SetPoint(N++,x,y,z);
   }
   for(Int_t i=0; i < 4000; i++){
      x = -25. + (generator->Rndm()*50.);
      y = -3.1 + (generator->Rndm()*6.2);
      z = 0;
      pm3h->SetPoint(N++,x,y,z);
   }
   for(Int_t i=0; i < 4000; i++){
      x = -25. + (generator->Rndm()*50.);
      y = 0;
      z = -60. + (generator->Rndm()*120.);
      pm3h->SetPoint(N++,x,y,z);
   }

   TVector3 v_Rot, v_Rot2, v_FARO;
   Double_t x_mean, y_mean, z_mean;
   Double_t x_meanT, y_meanT, z_meanT;
   Double_t x_G, y_G, z_G, r_G;
   Double_t x_dzn, y_dzn, z_dzn, phi_dzn, r_dzn;
   Double_t x_wid, y_wid, z_wid;
   Double_t x_max, x_min;
   Double_t y_max, y_min;
   Double_t z_max, z_min;
   Double_t x_mean_trans[256], y_mean_trans[256], z_mean_trans[256];
   TFile *inputfile = new TFile("/Users/taku/Downloads/survey_data/kiridashi.root","read");
   TFile *outputfile = new TFile("/Users/taku/Downloads/survey_data/Initial_value_DS.root","recreate");
   pm3->Write();
   pm3h->Write();
   TTree *tree = new TTree("tree","tree");
   tree->Branch("x_Local", &x_mean, "x_Local/D");
   tree->Branch("y_Local", &y_mean, "y_Local/D");
   tree->Branch("z_Local", &z_mean, "z_Local/D");
   tree->Branch("x_Global", &x_meanT, "x_Global/D");
   tree->Branch("y_Global", &y_meanT, "y_Global/D");
   tree->Branch("z_Global", &z_meanT, "z_Global/D");
   tree->Branch("x_Geo", &x_G, "x_Geo/D");
   tree->Branch("y_Geo", &y_G, "y_Geo/D");
   tree->Branch("z_Geo", &z_G, "z_Geo/D");
   tree->Branch("r_Geo", &r_G, "r_Geo/D");
   tree->Branch("x_dzn", &x_dzn, "x_dzn/D");
   tree->Branch("y_dzn", &y_dzn, "y_dzn/D");
   tree->Branch("z_dzn", &z_dzn, "z_dzn/D");
   tree->Branch("r_dzn", &r_dzn, "r_dzn/D");
   tree->Branch("phi_dzn", &phi_dzn, "phi_dzn/D");
   tree->Branch("x_wid", &x_wid, "x_wid/D");
   tree->Branch("y_wid", &y_wid, "y_wid/D");
   tree->Branch("z_wid", &z_wid, "z_wid/D");
   TTree *tree2 = new TTree("tree2","tree2");
   tree2->Branch("x_Local", &x_mean, "x_Local/D");
   tree2->Branch("y_Local", &y_mean, "y_Local/D");
   tree2->Branch("z_Local", &z_mean, "z_Local/D");
   tree2->Branch("x_Global", &x_meanT, "x_Global/D");
   tree2->Branch("y_Global", &y_meanT, "y_Global/D");
   tree2->Branch("z_Global", &z_meanT, "z_Global/D");
   tree2->Branch("x_Geo", &x_G, "x_Geo/D");
   tree2->Branch("y_Geo", &y_G, "y_Geo/D");
   tree2->Branch("z_Geo", &z_G, "z_Geo/D");
   tree2->Branch("r_Geo", &r_G, "r_Geo/D");
   tree2->Branch("x_dzn", &x_dzn, "x_dzn/D");
   tree2->Branch("y_dzn", &y_dzn, "y_dzn/D");
   tree2->Branch("z_dzn", &z_dzn, "z_dzn/D");
   tree2->Branch("r_dzn", &r_dzn, "r_dzn/D");
   tree2->Branch("phi_dzn", &phi_dzn, "phi_dzn/D");
   tree2->Branch("x_wid", &x_wid, "x_wid/D");
   tree2->Branch("y_wid", &y_wid, "y_wid/D");
   tree2->Branch("z_wid", &z_wid, "z_wid/D");
   TPolyMarker3D *pm3top[256], *pm3center[256];
   Int_t Niter = 0;
   TPolyMarker3D *pm3top_iter[256];
   TH1D *histX[256], *histY[256];
   TTree *transtree[256], *tree_copy[256];

   for(Int_t iPixel = 0; iPixel < 256; iPixel++){
      //if(iPixel%20) continue;
      Int_t USPixel = 511 - iPixel;
      cout << USPixel << endl;

      tree_copy[iPixel] = ((TTree*)inputfile->Get(Form("t3_%d",iPixel)))->CloneTree();
      if(!tree_copy[iPixel]){cout << "input tree cannot be read." << endl;}
      Int_t nEntry = tree_copy[iPixel]->GetEntries();
      tree_copy[iPixel]->SetBranchAddress("x",&x);
      tree_copy[iPixel]->SetBranchAddress("y",&y);
      tree_copy[iPixel]->SetBranchAddress("z",&z);
      tree_copy[iPixel]->SetBranchStatus("z2",0);
      tree_copy[iPixel]->SetBranchStatus("phi",0);
      tree_copy[iPixel]->SetBranchStatus("phi2",0);
      
      N = 0;
      histX[iPixel] = new TH1D(Form("histX_%d",USPixel),Form("histX_%d",USPixel),100,-35,-15);
      histY[iPixel] = new TH1D(Form("histY_%d",USPixel),Form("histY_%d",USPixel),100,-10,10);
      transtree[iPixel] = new TTree(Form("trans_%d",USPixel),Form("trans_%d",USPixel));
      transtree[iPixel]->Branch("x",&x,"x/D");
      transtree[iPixel]->Branch("y",&y,"y/D");
      transtree[iPixel]->Branch("z",&z,"z/D");
      pm3top[iPixel] = new TPolyMarker3D();
      pm3top[iPixel]->SetName(Form("pm3top_%d",USPixel));
      pm3top[iPixel]->SetMarkerStyle(8);
      pm3top[iPixel]->SetMarkerSize(0.1);
      pm3top[iPixel]->SetMarkerColor(kBlack);
   
      x_max = -40.;
      x_min = 40.;
      y_max = -10.;
      y_min = 10.;
      z_max = -75.;
      z_min = 75.;
      for(Int_t iEntry=0; iEntry < nEntry; iEntry++){
         tree_copy[iPixel]->GetEntry(iEntry);
         v_FARO = TVector3(x,y,z);
         v_Rot = v_FARO - v_par;
         v_Rot.RotateX(par[3]*TMath::DegToRad());
         v_Rot.RotateY(par[4]*TMath::DegToRad());
         v_Rot.RotateZ(par[5]*TMath::DegToRad());
         v_Rot2 = v_Rot - v_geo[USPixel];
         v_Rot2.RotateZ(-psi_geo[USPixel]);
         v_Rot2.RotateX(-theta_geo[USPixel]);
         v_Rot2.RotateZ(-phi_geo[USPixel]);
         x = v_Rot2.X();
         y = v_Rot2.Y();
         z = v_Rot2.Z();
         transtree[iPixel]->Fill();
         if(x < -19 && -22 < x){
            if((USPixel >= 288 && USPixel <= 302) && USPixel%2 == 0) continue;
            else if((USPixel >= 384 && USPixel <= 399) && USPixel%2 != 0) continue;
            else if(USPixel >= 304 && USPixel <= 383) continue;
            pm3top[iPixel]->SetPoint(N++,x,y,z);
            histX[iPixel]->Fill(x);
            histY[iPixel]->Fill(y);
         }
         if(x < -24 && -27 < x){
            if((USPixel >= 288 && USPixel <= 302) && USPixel%2 != 0) continue;
            else if((USPixel >= 384 && USPixel <= 399) && USPixel%2 == 0) continue;
            else if(USPixel <= 287 && USPixel >= 400) continue;
            pm3top[iPixel]->SetPoint(N++,x,y,z);
            histX[iPixel]->Fill(x);
            histY[iPixel]->Fill(y);
         }
      }
      transtree[iPixel]->Write();
      if(N!= 0){
         Niter = 0;
         pm3top_iter[iPixel] = new TPolyMarker3D();
         pm3top_iter[iPixel]->SetName(Form("pm3top_iter_%d",USPixel));
         pm3top_iter[iPixel]->SetMarkerStyle(8);
         pm3top_iter[iPixel]->SetMarkerSize(0.15);
         pm3top_iter[iPixel]->SetMarkerColor(kRed);
         for(Int_t i=0; i<N; i++){
            pm3top[iPixel]->GetPoint(i,x,y,z);
            if(fabs(x - histX[iPixel]->GetMean()) > histX[iPixel]->GetRMS()) continue;
            if(fabs(y - histY[iPixel]->GetMean()) > histY[iPixel]->GetRMS()*2) continue;
            if(fabs(z) > 75.) continue;
            pm3top_iter[iPixel]->SetPoint(Niter++,x,y,z);
            if(x > x_max){x_max = x;}
            if(x < x_min){x_min = x;}
            if(y > y_max){y_max = y;}
            if(y < y_min){y_min = y;}
            if(z > z_max){z_max = z;}
            if(z < z_min){z_min = z;}
         }
         pm3center[iPixel] = new TPolyMarker3D();
         pm3center[iPixel]->SetName(Form("pm3center_%d",USPixel));
         pm3center[iPixel]->SetMarkerStyle(8);
         pm3center[iPixel]->SetMarkerSize(0.15);
         pm3center[iPixel]->SetMarkerColor(kGreen+2);
         x_mean = (x_max + x_min)/2 + 20.;
         if((USPixel >= 288 && USPixel <= 302) && USPixel%2 == 0) x_mean += 5.;
         else if((USPixel >= 384 && USPixel <= 399) && USPixel%2 != 0) x_mean +=5.;
         else if(USPixel >= 304 && USPixel <= 383) x_mean +=5.;
         y_mean = (y_max + y_min)/2;
         z_mean = (z_max + z_min)/2;
         x_dzn = v_geo[USPixel].X();
         y_dzn = v_geo[USPixel].Y();
         z_dzn = v_geo[USPixel].Z();
         r_dzn = sqrt(x_dzn * x_dzn + y_dzn * y_dzn);
         phi_dzn = TMath::ATan(y_dzn/x_dzn) * TMath::RadToDeg();
         if(phi_dzn > 0.){phi_dzn -= 180.;}
         x_wid = x_max - x_min;
         y_wid = y_max - y_min;
         z_wid = z_max - z_min;
         v_Rot = TVector3(x_mean, y_mean, z_mean);
         v_Rot.RotateZ(phi_geo[USPixel]);
         v_Rot.RotateX(theta_geo[USPixel]);
         v_Rot.RotateZ(psi_geo[USPixel]);
         x_meanT = v_Rot.X();
         y_meanT = v_Rot.Y();
         z_meanT = v_Rot.Z();
         x_G = x_dzn + x_meanT;
         y_G = y_dzn + y_meanT;
         z_G = z_dzn + z_meanT;
         r_G = sqrt(x_G * x_G + y_G * y_G);
         tree->Fill();
         x_mean_trans[255-iPixel] = v_Rot.X();
         y_mean_trans[255-iPixel] = v_Rot.Y();
         z_mean_trans[255-iPixel] = v_Rot.Z();
         std::cout << x_mean << ", " << y_mean << ", " << z_mean <<
             ", "  << x_wid << ", " << y_wid << ", " << z_wid <<  std::endl;
         for(Int_t iPoint=0; iPoint<3000; iPoint++){
            x = -20. + (generator->Rndm()*40.) + x_mean;
            y = y_mean;
            z = -60. + (generator->Rndm()*120.) + z_mean;
            if((USPixel >= 288 && USPixel <= 302) && USPixel%2 == 0) x *= 1.25;
            else if((USPixel >= 384 && USPixel <= 399) && USPixel%2 != 0) x *= 1.25;
            else if(USPixel >= 304 && USPixel <= 383) x *= 1.25;
            pm3center[iPixel]->SetPoint(iPoint,x,y,z);
         }
         for(Int_t iPoint=3000; iPoint<6000; iPoint++){
            x = -20. + (generator->Rndm()*40.) + x_mean;
            y = -3.1 + (generator->Rndm()*6.2) + y_mean;
            z = z_mean;
            if((USPixel >= 288 && USPixel <= 302) && USPixel%2 == 0) x *= 1.25;
            else if((USPixel >= 384 && USPixel <= 399) && USPixel%2 != 0) x *= 1.25;
            else if(USPixel >= 304 && USPixel <= 383) x *= 1.25;
            pm3center[iPixel]->SetPoint(iPoint,x,y,z);
         }
         for(Int_t iPoint=6000; iPoint<9000; iPoint++){
            x = x_mean;
            y = -3.1 + (generator->Rndm()*6.2) + y_mean;
            z = -60. + (generator->Rndm()*120.) + z_mean;
            pm3center[iPixel]->SetPoint(iPoint,x,y,z);
         }
         pm3top[iPixel]->Write();
         pm3top_iter[iPixel]->Write();
         pm3center[iPixel]->Write();
      }
   }
   tree->Write();
   vector<Int_t> ex_id;
   for(Int_t i=0; i<256; i++){
      tree->GetEntry(i);
      if(fabs(y_wid - 7.1409369) < 2*0.68583155 && fabs(z_wid - 126.69637) < 2*3.5845383){
         if(i==47 || i==63 || i==95 || i==99 || i==126 || i==127 || i==143 || i==158 || i==159 || i==169){
            ex_id.push_back(i);
            x_mean = 0.;
            y_mean = 0.;
            z_mean = 0.;
            x_meanT = 0.;
            y_meanT = 0.;
            z_meanT = 0.;
            x_G = x_dzn;
            y_G = y_dzn;
            z_G = z_dzn;
            r_G = r_dzn;
            tree2->Fill();
         }else{
            tree2->Fill();
         }
      }else{
         ex_id.push_back(i);
         x_mean = 0.;
         y_mean = 0.;
         z_mean = 0.;
         x_meanT = 0.;
         y_meanT = 0.;
         z_meanT = 0.;
         x_G = x_dzn;
         y_G = y_dzn;
         z_G = z_dzn;
         r_G = r_dzn;
         tree2->Fill();
      }
   }
   tree2->Write();
   outputfile->Close();

   /*
   Int_t ex_flag = 0;
   Int_t ex_size = ex_id.size();
   std::ofstream geometry2("outputfiles/Geometry2.csv");
   geometry2 << "id \t idx \t XYZPosition_0 \t XYZPosition_1 \t XYZPosition_2 \t XPhi \t XTheta \t XPsi \t Order \t GeometricalHitOrder" << std::endl;
   for(Int_t i=0; i<512; i++){
      if(i < 256){
         geometry2 << id_geo[i]+1 << "\t" << idx_geo[i] << "\t" << v_geo[i].X()/10. << "\t" << v_geo[i].Y()/10. <<
            "\t" << v_geo[i].Z()/10. << "\t" << phi_geo[i] << "\t" << theta_geo[i] << "\t" << psi_geo[i] <<
            "\t" << order_geo[i] << "\t" << hitorder_geo[i] << std::endl;
      }
      if(i >= 256){
         for(Int_t j=0; j < ex_size; j++){
            if(511 - i == ex_id[j]) ex_flag = 1;
         }
         if(ex_flag){
            geometry2 << id_geo[i]+1 << "\t" << idx_geo[i] << "\t" << v_geo[i].X()/10. << "\t" << v_geo[i].Y()/10. <<
               "\t" << v_geo[i].Z()/10. << "\t" << phi_geo[i] << "\t" << theta_geo[i] << "\t" << psi_geo[i] <<
               "\t" << order_geo[i] << "\t" << hitorder_geo[i] << std::endl;
         }else{
            geometry2 << id_geo[i]+1 << "\t" << idx_geo[i] << "\t" << (v_geo[i].X() - x_mean_trans[i-256])/10. <<
            "\t" << (v_geo[i].Y() - y_mean_trans[i-256])/10. << "\t" << (v_geo[i].Z() - z_mean_trans[i-256])/10. <<
            "\t" << phi_geo[i] << "\t" << theta_geo[i] << "\t" << psi_geo[i] <<
            "\t" << order_geo[i] << "\t" << hitorder_geo[i] << std::endl;
         }
         ex_flag = 0;
      }
   }
   geometry2.close();

   Double_t Rndm[1536];
   std::ofstream PixelPosition("outputfiles/PixelPosition_Rndm5.dat");
   for(Int_t i=0; i<512; i++){
      if(i < 256){
         Rndm[3*i] = (generator->Rndm()-0.5)*10;
         Rndm[3*i+1] = (generator->Rndm()-0.5)*10;
         Rndm[3*i+2] = (generator->Rndm()-0.5)*10;
         x = v_geo[i].X() + Rndm[3*i];
         y = v_geo[i].Y() + Rndm[3*i+1];
         z = v_geo[i].Z() + Rndm[3*i+2];
         phi = TMath::ATan(y/x) * TMath::RadToDeg();
         if(phi > 0){phi -= 180.;}
         z /= 10.;
         PixelPosition << idx_geo[i] << "\t";
         PixelPosition << std::fixed << std::setw(6) << z << "\t" << phi << std::endl;
      }
      if(i >= 256){
         Rndm[3*i] = (generator->Rndm()-0.5)*10;
         Rndm[3*i+1] = (generator->Rndm()-0.5)*10;
         Rndm[3*i+2] = (generator->Rndm()-0.5)*10;
         x = v_geo[i].X() + Rndm[3*i];
         y = v_geo[i].Y() + Rndm[3*i+1];
         z = v_geo[i].Z() + Rndm[3*i+2];
         phi = TMath::ATan(y/x) * TMath::RadToDeg();
         if(phi > 0){phi -= 180.;}
         z /= 10.;
         PixelPosition << idx_geo[i] << "\t";
         PixelPosition << std::fixed << std::setw(6) << z << "\t" << phi << std::endl;
      }
   }
   PixelPosition.close();
   
   std::ofstream Rndms("outputfiles/Rndm5_dat");
   for(Int_t i=0; i<1536; i++){
      Rndms << Rndm[3*i] << "\t" << Rndm[3*i+1] << "\t" << Rndm[3*i+2] << std::endl;
   }
   Rndms.close();
   */
}

