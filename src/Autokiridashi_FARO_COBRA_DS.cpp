#include <TROOT.h>
#include <TMath.h>
#include <TVector3.h>
#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TGraph2D.h>
#include <TPolyMarker3D.h>
#include <TRandom.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;


void Autokiridashi_FARO_COBRA_DS()
{
   vector<Double_t> x_geo, y_geo, z_geo;
   vector<Double_t> phi_geo, theta_geo, psi_geo;
   vector<Double_t> par, par_err;
   vector<Int_t> id_geo,idx_geo,order_geo,hitorder_geo;
   Double_t x,y,z,phi,theta,psi;
   
   Int_t id,idx,order,hitorder;
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
      x_geo.push_back(x);
      y_geo.push_back(y);
      z_geo.push_back(z);
      phi_geo.push_back(phi);
      theta_geo.push_back(theta);
      psi_geo.push_back(psi);
      id_geo.push_back(id);
      idx_geo.push_back(idx);
      order_geo.push_back(order);
      hitorder_geo.push_back(hitorder);
   }
   ifs.close();

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
   pm3->SetName("pm3_origin");
   pm3->SetMarkerStyle(8);
   pm3->SetMarkerSize(0.1);
   pm3->SetMarkerColor(kBlue);
   Int_t N = 0;
   for(Int_t i=0; i < 3000; i++){
      x = 0;
      y = 0;
      z = -60. + (generator->Rndm()*120.);
      pm3->SetPoint(N++,x,y,z);
   }
   for(Int_t i=0; i < 2000; i++){
      x = 0;
      y = -6. + (generator->Rndm()*12.);
      z = 0;
      pm3->SetPoint(N++,x,y,z);
   }
   for(Int_t i=0; i < 3000; i++){
      x = -20. + (generator->Rndm()*40.);
      y = 0;
      z = 0;
      pm3->SetPoint(N++,x,y,z);
   }

   
   TVector3 v;
   Double_t a,b,c,xt,yt,zt;
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

   TFile *fout = new TFile("/Users/taku/Downloads/survey_data/Autokiridashi_FARO_COBRA_DS.root","recreate");
   TPolyMarker3D *pm3top[256], *pm3center[256];
   TPolyMarker3D *pm3top_iter[256];
   TTree *td = new TTree("data","data");
   td->Branch("x",&x,"x/D"); 
   td->Branch("y",&y,"y/D");
   td->Branch("z",&z,"z/D");
   ifs.open("/Users/taku/Downloads/survey_data/20190320_MEG_II_Timing_Counter_Downstream_Final.asc");
   while(ifs >> x >> y >> z >> a >> b >> c) td->Fill();
   ifs.close();
   td->Write();
   pm3->Write();

   fout->cd();
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

   Int_t Niter = 0;
   TH1D *histX[256], *histY[256];
   Int_t nEntry = td->GetEntries();

   for(Int_t iPixel = 0; iPixel < 256; iPixel++){
      TGraph2D *g2d = new TGraph2D();
      N = 0;
      cout << iPixel << endl;
      Int_t DSPixel = iPixel;
      TTree *t3 = new TTree(Form("t3_%d",iPixel),Form("t3_%d",iPixel));
      t3->Branch("x",&x,"x/D");
      t3->Branch("y",&y,"y/D");
      t3->Branch("z",&z,"z/D");
      t3->Branch("xt",&xt,"xt/D");
      t3->Branch("yt",&yt,"yt/D");
      t3->Branch("zt",&zt,"zt/D");

      histX[iPixel] = new TH1D(Form("histX_%d",DSPixel),Form("histX_%d",DSPixel),100,-35,-15);
      histY[iPixel] = new TH1D(Form("histY_%d",DSPixel),Form("histY_%d",DSPixel),100,-10,10);

      pm3top[iPixel] = new TPolyMarker3D();
      pm3top[iPixel]->SetName(Form("pm3top_%d",DSPixel));
      pm3top[iPixel]->SetMarkerStyle(8);
      pm3top[iPixel]->SetMarkerSize(0.1);
      pm3top[iPixel]->SetMarkerColor(kBlack);
      
      for(Int_t iEntry = 0; iEntry < nEntry; iEntry++){
         td->GetEntry(iEntry);
         v_Rot = TVector3(x,y,z);
         TVector3 v_geo(x_geo[DSPixel],y_geo[DSPixel],z_geo[DSPixel]);
         v_Rot = v_Rot - v_par;
         v_Rot.RotateX(par[3]*TMath::DegToRad());
         v_Rot.RotateY(par[4]*TMath::DegToRad());
         v_Rot.RotateZ(par[5]*TMath::DegToRad());
         v_Rot = v_Rot - v_geo;
         v_Rot.RotateZ(-psi_geo[DSPixel]);
         v_Rot.RotateX(-theta_geo[DSPixel]);
         v_Rot.RotateZ(-phi_geo[DSPixel]);
         xt = v_Rot.X();
         yt = v_Rot.Y();
         zt = v_Rot.Z();
         if(fabs(xt) < 27. && fabs(yt) < 10. && fabs(zt) < 80.){
           t3->Fill();
           if(xt < -19 && -22 < xt){
              if((DSPixel >= 32 && DSPixel <= 46) && DSPixel%2 == 0) continue;
              else if((DSPixel >= 128 && DSPixel <= 143) && DSPixel%2 != 0) continue;
              else if(DSPixel >= 48 && DSPixel <= 127) continue;
              pm3top[iPixel]->SetPoint(N++,xt,yt,zt);
              histX[iPixel]->Fill(xt);
              histY[iPixel]->Fill(yt);
           }
           if(xt < -24 && -27 < xt){
              if((DSPixel >= 32 && DSPixel <= 46) && DSPixel%2 != 0) continue;
              else if((DSPixel >= 128 && DSPixel <= 143) && DSPixel%2 == 0) continue;
              else if(DSPixel <= 31 && DSPixel >= 144) continue;
              pm3top[iPixel]->SetPoint(N++,xt,yt,zt);
              histX[iPixel]->Fill(xt);
              histY[iPixel]->Fill(yt);
           }
         }
      }
      
      t3->Write();
      
      x_max = -40.;
      x_min = 40.;
      y_max = -10.;
      y_min = 10.;
      z_max = -75.;
      z_min = 75.;
      
      if(N!= 0){
         Niter = 0;
         pm3top_iter[iPixel] = new TPolyMarker3D();
         pm3top_iter[iPixel]->SetName(Form("pm3top_iter_%d",DSPixel));
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
         pm3center[iPixel]->SetName(Form("pm3center_%d",DSPixel));
         pm3center[iPixel]->SetMarkerStyle(8);
         pm3center[iPixel]->SetMarkerSize(0.15);
         pm3center[iPixel]->SetMarkerColor(kGreen+2);

         x_mean = (x_max + x_min)/2 + 20.;
         if((DSPixel >= 32 && DSPixel <= 46) && DSPixel%2 == 0) x_mean += 5.;
         else if((DSPixel >= 128 && DSPixel <= 143) && DSPixel%2 != 0) x_mean +=5.;
         else if(DSPixel >= 48 && DSPixel <= 127) x_mean +=5.;
         y_mean = (y_max + y_min)/2;
         z_mean = (z_max + z_min)/2;
         x_dzn = x_geo[DSPixel];
         y_dzn = y_geo[DSPixel];
         z_dzn = z_geo[DSPixel];
         r_dzn = sqrt(x_dzn * x_dzn + y_dzn * y_dzn);
         phi_dzn = TMath::ATan(y_dzn/x_dzn) * TMath::RadToDeg();
         if(phi_dzn > 0.){phi_dzn -= 180.;}
         x_wid = x_max - x_min;
         y_wid = y_max - y_min;
         z_wid = z_max - z_min;
         v_Rot = TVector3(x_mean, y_mean, z_mean);
         v_Rot.RotateZ(phi_geo[DSPixel]);
         v_Rot.RotateX(theta_geo[DSPixel]);
         v_Rot.RotateZ(psi_geo[DSPixel]);
         x_meanT = v_Rot.X();
         y_meanT = v_Rot.Y();
         z_meanT = v_Rot.Z();
         x_G = x_dzn + x_meanT;
         y_G = y_dzn + y_meanT;
         z_G = z_dzn + z_meanT;
         r_G = sqrt(x_G * x_G + y_G * y_G);
         tree->Fill();
         x_mean_trans[iPixel] = v_Rot.X();
         y_mean_trans[iPixel] = v_Rot.Y();
         z_mean_trans[iPixel] = v_Rot.Z();
         std::cout << x_mean << ", " << y_mean << ", " << z_mean <<
             ", "  << x_wid << ", " << y_wid << ", " << z_wid <<  std::endl;
         for(Int_t iPoint=0; iPoint<3000; iPoint++){
            x = -20. + (generator->Rndm()*40.) + x_mean;
            y = y_mean;
            z = -60. + (generator->Rndm()*120.) + z_mean;
            if((DSPixel >= 32 && DSPixel <= 46) && DSPixel%2 == 0) x *= 1.25;
            else if((DSPixel >= 128 && DSPixel <= 143) && DSPixel%2 != 0) x *= 1.25;
            else if(DSPixel >= 48 && DSPixel <= 127) x *= 1.25;
            pm3center[iPixel]->SetPoint(iPoint,x,y,z);
         }
         for(Int_t iPoint=3000; iPoint<6000; iPoint++){
            x = -20. + (generator->Rndm()*40.) + x_mean;
            y = -3.1 + (generator->Rndm()*6.2) + y_mean;
            z = z_mean;
            if((DSPixel >= 32 && DSPixel <= 46) && DSPixel%2 == 0) x *= 1.25;
            else if((DSPixel >= 128 && DSPixel <= 143) && DSPixel%2 != 0) x *= 1.25;
            else if(DSPixel >= 48 && DSPixel <= 127) x *= 1.25;
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
      if(fabs(y_wid - 6.5018153) < 2*1.0207279 && fabs(z_wid - 122.69740) < 2*4.0111153){
         tree2->Fill();
      }else{
         ex_id.push_back(i);
         x_mean = 0.;
         y_mean = 0.;
         z_mean = 0.;
         x_meanT = 0.;
         y_meanT = 0.;
         z_meanT = 0.;
         x_wid = 0.;
         y_wid = 0.;
         z_wid = 0.;
         x_G = x_dzn;
         y_G = y_dzn;
         z_G = z_dzn;
         r_G = r_dzn;
         tree2->Fill();
      }
   }
   tree2->Write();
   fout->Close();
   
   Int_t ex_flag = 0;
   Int_t ex_size = ex_id.size();
   std::ofstream geometry2("outputfiles/Geometry2_DS.csv");
   geometry2 << "id \t idx \t XYZPosition_0 \t XYZPosition_1 \t XYZPosition_2 \t XPhi \t XTheta \t XPsi \t Order \t GeometricalHitOrder" << std::endl;
   for(Int_t i=0; i<512; i++){
      if(i >= 256){
         geometry2 << id_geo[i]+1 << "\t" << idx_geo[i] << "\t" << x_geo[i]/10. << "\t" << y_geo[i]/10. <<
            "\t" << z_geo[i]/10. << "\t" << phi_geo[i] << "\t" << theta_geo[i] << "\t" << psi_geo[i] <<
            "\t" << order_geo[i] << "\t" << hitorder_geo[i] << std::endl;
      }
      if(i < 256){
         for(Int_t j=0; j < ex_size; j++){
            if(i == ex_id[j]) ex_flag = 1;
         }
         if(ex_flag){
            geometry2 << id_geo[i]+1 << "\t" << idx_geo[i] << "\t" << x_geo[i]/10. << "\t" << y_geo[i]/10. <<
               "\t" << z_geo[i]/10. << "\t" << phi_geo[i] << "\t" << theta_geo[i] << "\t" << psi_geo[i] <<
               "\t" << order_geo[i] << "\t" << hitorder_geo[i] << std::endl;
         }else{
            geometry2 << id_geo[i]+1 << "\t" << idx_geo[i] << "\t" << (x_geo[i] - x_mean_trans[i])/10. <<
            "\t" << (y_geo[i] - y_mean_trans[i])/10. << "\t" << (z_geo[i] - z_mean_trans[i])/10. <<
            "\t" << phi_geo[i] << "\t" << theta_geo[i] << "\t" << psi_geo[i] <<
            "\t" << order_geo[i] << "\t" << hitorder_geo[i] << std::endl;
         }
         ex_flag = 0;
      }
   }
   geometry2.close();

   ex_flag = 0;
   std::ofstream PixelPosition("outputfiles/PixelPosition_2018DSsurvey.dat");
   for(Int_t i=0; i<512; i++){
      for(Int_t j=0; j < ex_size; j++){
         if(i == ex_id[j]) ex_flag = 1;
      }
      if(i < 256 && ex_flag == 0){
         x = x_geo[i] - x_mean_trans[i];
         y = y_geo[i] - y_mean_trans[i];
         z = z_geo[i] - z_mean_trans[i];
         phi = TMath::ATan(y/x) * TMath::RadToDeg();
         if(phi > 0){phi -= 180.;}
         z /= 10.;
         PixelPosition << idx_geo[i] << "\t";
         PixelPosition << std::fixed << std::setw(6) << z << "\t" << phi << std::endl;
      }
      else if(i >= 256 || ex_flag == 1){
         x = x_geo[i];
         y = y_geo[i];
         z = z_geo[i];
         phi = TMath::ATan(y/x) * TMath::RadToDeg();
         if(phi > 0){phi -= 180.;}
         z /= 10.;
         PixelPosition << idx_geo[i] << "\t";
         PixelPosition << std::fixed << std::setw(6) << z << "\t" << phi << std::endl;
      }
      ex_flag = 0;
   }
   PixelPosition.close();

}

