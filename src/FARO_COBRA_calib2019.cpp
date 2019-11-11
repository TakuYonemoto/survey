#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TMath.h>
#include <TMinuit.h>

using namespace TMath;

void LossFunc(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t flag);
void check(Double_t *par);

TVector3 v3_FARO = TVector3(-2.7837, -0.0089, 0.0279);
TVector3 v9_FARO = TVector3(457.4944, 351.0098, 152.2236);
TVector3 v10_FARO = TVector3(391.2457, 362.3315, 148.6683);
TVector3 v11_FARO = TVector3(139.8866, 290.4434, 148.947);

TVector3 v7_COBRA = TVector3(-358.801, -75.263, -1147.147);
TVector3 v9_COBRA = TVector3(158.605, -338.287, -1294.136);
TVector3 v10_COBRA = TVector3(95.281, -361.259, -1291.291);
TVector3 v11_COBRA = TVector3(-164.731, -335.467, -1294.349);

Int_t npar = 6;

void FARO_COBRA_calib2019(){

   /*
// Geometrical Method    
   TVector3 v3_F, v9_F, v11_F;
   TVector3 v7_C, v9_C, v11_C;

   v3_F = v3_FARO - v10_FARO;
   v7_C = v7_COBRA - v10_COBRA;
   v9_F = v9_FARO - v10_FARO;
   v11_F = v11_FARO - v10_FARO;
   v9_C = v9_COBRA - v10_COBRA;
   v11_C = v11_COBRA - v10_COBRA;

   //-----rotate-----
   TVector3 v_Rot, v_Tar, v_Ax, v_Ax2;
   
   v_Rot = v9_F * (1./v9_F.Mag());
   v_Tar = v9_C * (1./v9_C.Mag());
   v_Ax2 = v_Tar;
   Double_t theta9 = ACos(v_Rot * v_Tar);
   v_Ax = v_Rot.Cross(v_Tar);
   v_Rot = v11_F * (1./v11_F.Mag());
   v_Tar = v11_C * (1./v11_C.Mag());
   v_Rot.Rotate(theta9, v_Ax);
   Double_t theta11_about_v9 = 2*ATan(((v_Rot - v_Tar).Mag())/((v_Ax2.Cross(v_Rot + v_Tar)).Mag()));

   //----check-----
   TVector3 v_Geo;
   Double_t Loss = 0.;
   
   v_Rot = v3_F;
   v_Rot.Rotate(theta9, v_Ax);
   v_Rot.Rotate(theta11_about_v9, v_Ax2);
   v_Geo = v_Rot - v7_C;
   std::cout << " 3 : " << v_Geo.X() << " , " << v_Geo.Y() << " , " << v_Geo.Z() << std::endl;
   Loss += v_Geo.Mag2();
   v_Rot = v9_F;
   v_Rot.Rotate(theta9, v_Ax);
   v_Rot.Rotate(theta11_about_v9, v_Ax2);
   v_Geo = v_Rot - v9_C;
   std::cout << " 9 : " << v_Geo.X() << " , " << v_Geo.Y() << " , " << v_Geo.Z() << std::endl;
   Loss += v_Geo.Mag2();
   v_Rot = v11_F;
   v_Rot.Rotate(theta9, v_Ax);
   v_Rot.Rotate(theta11_about_v9, v_Ax2);
   v_Geo = v_Rot - v11_C;
   std::cout << " 11 : " << v_Geo.X() << " , " << v_Geo.Y() << " , " << v_Geo.Z() << std::endl;
   Loss += v_Geo.Mag2();
   std::cout << "Loss = " << Loss << std::endl;
   */

// Fitting Method
   TMinuit *minuit = new TMinuit(npar);
   minuit->SetFCN(LossFunc);
   Int_t ierflg = 0;
   const char *parname[] = {"x_parallel","y_parallel","z_parallel","theta_x","theta_y","theta_z"};
   Double_t vstart[] = {376.,-9.,-1143.,359,180,190};
   
   for(Int_t i=0; i < npar; i++){
      minuit->mnparm(i, parname[i],vstart[i], 0.1, 0, 0, ierflg);
   }

   Double_t arglist[10];
   arglist[0] = TMath::ChisquareQuantile(0.682689492,npar);
   minuit->mnexcm("SET ERR", arglist, 1, ierflg);

   arglist[0] = 10000;
   arglist[1] = 1;
   minuit->mnexcm("MINOS", arglist, 2, ierflg);
   Double_t fedm, errdef, fmin;
   Int_t npari, nparx, istat;
   minuit->mnstat(fmin,fedm,errdef,npari,nparx,istat);

   //Double_t redchi = fmin/(4*3 -6);
   std::cout << "LossFunc = " << fmin << std::endl;
   Double_t par_result[npar], par_err[npar];
   for(Int_t i=0; i<npar; i++){
      Double_t current, error;
      minuit->GetParameter(i, current, error);
      par_result[i] = current;
      par_err[i] = error;
      if(i==0) std::cout << "parname : value[mm] : error[mm] " << std::endl;
      else if(i==3) std::cout << "parname : value[degree] : error[degree] " << std::endl;
      std::cout << parname[i] << " : " << current << " : " << error << std::endl;
   }

   std::cout << "=======================================" << std::endl;
   check(par_result);
   std::ofstream output_file("outputfiles/FARO_COBRA_calib2019.csv");
   for(Int_t i=0; i < npar; i++){
      output_file << par_result[i] << "," << par_err[i] << std::endl;
   }
   output_file.close();
}


void LossFunc(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t flag)
{
   Double_t L = 0.;

   TVector3 v_Rot,v_Res;
   TVector3 v_Para(par[0], par[1], par[2]);
   v_Rot = v3_FARO - v_Para;
   v_Rot.RotateX(par[3]*TMath::DegToRad());
   v_Rot.RotateY(par[4]*TMath::DegToRad());
   v_Rot.RotateZ(par[5]*TMath::DegToRad());
   v_Res = v7_COBRA - v_Rot;
   L += v_Res.Mag2();

   v_Rot = v9_FARO - v_Para;
   v_Rot.RotateX(par[3]*TMath::DegToRad());
   v_Rot.RotateY(par[4]*TMath::DegToRad());
   v_Rot.RotateZ(par[5]*TMath::DegToRad());
   v_Res = v9_COBRA - v_Rot;
   L += v_Res.Mag2();

   v_Rot = v10_FARO - v_Para;
   v_Rot.RotateX(par[3]*TMath::DegToRad());
   v_Rot.RotateY(par[4]*TMath::DegToRad());
   v_Rot.RotateZ(par[5]*TMath::DegToRad());
   v_Res = v10_COBRA - v_Rot;
   L += v_Res.Mag2();
   
   v_Rot = v11_FARO - v_Para;
   v_Rot.RotateX(par[3]*TMath::DegToRad());
   v_Rot.RotateY(par[4]*TMath::DegToRad());
   v_Rot.RotateZ(par[5]*TMath::DegToRad());
   v_Res = v11_COBRA - v_Rot;
   L += v_Res.Mag2();

   f = L;

}

void check(Double_t *par)
{
   Double_t L = 0;
   
   TVector3 v_Rot,v_Res;
   TVector3 v_Para(par[0], par[1], par[2]);
   v_Rot = v3_FARO - v_Para;
   v_Rot.RotateX(par[3]*TMath::DegToRad());
   v_Rot.RotateY(par[4]*TMath::DegToRad());
   v_Rot.RotateZ(par[5]*TMath::DegToRad());
   v_Res = v7_COBRA - v_Rot;
   std::cout << " 3 : " << v_Res.X() << ", " << v_Res.Y() << ", " << v_Res.Z() << std::endl;
   L += v_Res.Mag2();

   v_Rot = v9_FARO - v_Para;
   v_Rot.RotateX(par[3]*TMath::DegToRad());
   v_Rot.RotateY(par[4]*TMath::DegToRad());
   v_Rot.RotateZ(par[5]*TMath::DegToRad());
   v_Res = v9_COBRA - v_Rot;
   std::cout << " 9 : " << v_Res.X() << ", " << v_Res.Y() << ", " << v_Res.Z() << std::endl;
   L += v_Res.Mag2();

   v_Rot = v10_FARO - v_Para;
   v_Rot.RotateX(par[3]*TMath::DegToRad());
   v_Rot.RotateY(par[4]*TMath::DegToRad());
   v_Rot.RotateZ(par[5]*TMath::DegToRad());
   v_Res = v10_COBRA - v_Rot;
   std::cout << " 10 : " << v_Res.X() << ", " << v_Res.Y() << ", " << v_Res.Z() << std::endl;
   L += v_Res.Mag2();
   
   v_Rot = v11_FARO - v_Para;
   v_Rot.RotateX(par[3]*TMath::DegToRad());
   v_Rot.RotateY(par[4]*TMath::DegToRad());
   v_Rot.RotateZ(par[5]*TMath::DegToRad());
   v_Res = v11_COBRA - v_Rot;
   std::cout << " 11 : " << v_Res.X() << ", " << v_Res.Y() << ", " << v_Res.Z() << std::endl;
   L += v_Res.Mag2();

   std::cout << "LossFunc = " << L << endl;
 
   //1.5
   TVector3 old1_FARO = TVector3(751.7787, -23.3314, 79.8673);
   TVector3 old2_FARO = TVector3(720.1456, -23.0381, -44.4701);
   TVector3 old3_FARO = TVector3(-26.7815, -26.9501, 69.5445);
   TVector3 old4_FARO = TVector3(7.3234, -26.9042, -53.1685);
   //0.5
   TVector3 old1_FARO2 = TVector3(751.9625, -23.9373, 65.4922);
   TVector3 old2_FARO2 = TVector3(734.589, -23.7316, -44.8804);
   TVector3 old3_FARO2 = TVector3(-26.5333, -26.9665, 55.0457);
   TVector3 old4_FARO2 = TVector3(-7.1412, -26.8283, -53.4788);
   //Markers
   TVector3 old1_CAD = TVector3(380.17, 83.901, -1218.4);   //Marker A
   TVector3 old3_CAD = TVector3(-385.94, -51.184, -1218.4); //Marker B
   TVector3 old2_CAD = TVector3(348.51, 78.733, -1092);     //Marker C
   TVector3 old4_CAD = TVector3(-354.42, -45.213, -1092);   //Marker D
   
   std::cout<<"      ID       : x_diff, y_diff, z_diff ([mm] in COBRA coodinate system)"<<std::endl;
   v_Rot = old1_FARO - v_Para;
   v_Rot.RotateX(par[3]*TMath::DegToRad());
   v_Rot.RotateY(par[4]*TMath::DegToRad());
   v_Rot.RotateZ(par[5]*TMath::DegToRad());
   v_Res = v_Rot - old1_CAD;
   std::cout<<" old1_1.5_diff : "<< v_Res.X()<<", "<<v_Res.Y()<<", "<<v_Res.Z()<<std::endl;
   
   v_Rot = old2_FARO - v_Para;
   v_Rot.RotateX(par[3]*TMath::DegToRad());
   v_Rot.RotateY(par[4]*TMath::DegToRad());
   v_Rot.RotateZ(par[5]*TMath::DegToRad());
   v_Res = v_Rot - old2_CAD;
   std::cout<<" old2_1.5_diff : "<< v_Res.X()<<", "<<v_Res.Y()<<", "<<v_Res.Z()<<std::endl;

   v_Rot = old3_FARO - v_Para;
   v_Rot.RotateX(par[3]*TMath::DegToRad());
   v_Rot.RotateY(par[4]*TMath::DegToRad());
   v_Rot.RotateZ(par[5]*TMath::DegToRad());
   v_Res = v_Rot - old3_CAD;
   std::cout<<" old3_1.5_diff : "<< v_Res.X()<<", "<<v_Res.Y()<<", "<<v_Res.Z()<<std::endl;

   v_Rot = old4_FARO - v_Para;
   v_Rot.RotateX(par[3]*TMath::DegToRad());
   v_Rot.RotateY(par[4]*TMath::DegToRad());
   v_Rot.RotateZ(par[5]*TMath::DegToRad());
   v_Res = v_Rot - old4_CAD;
   std::cout<<" old4_1.5_diff : "<< v_Res.X()<<", "<<v_Res.Y()<<", "<<v_Res.Z()<<std::endl;

   std::cout<<"***************************************************"<<std::endl;
   v_Rot = old1_FARO2 - v_Para;
   v_Rot.RotateX(par[3]*TMath::DegToRad());
   v_Rot.RotateY(par[4]*TMath::DegToRad());
   v_Rot.RotateZ(par[5]*TMath::DegToRad());
   v_Res = v_Rot - old1_CAD;
   std::cout<<" old1_0.5_diff : "<< v_Res.X()<<", "<<v_Res.Y()<<", "<<v_Res.Z()<<std::endl;
   v_Rot = old2_FARO2 - v_Para;
   v_Rot.RotateX(par[3]*TMath::DegToRad());
   v_Rot.RotateY(par[4]*TMath::DegToRad());
   v_Rot.RotateZ(par[5]*TMath::DegToRad());
   v_Res = v_Rot - old2_CAD;
   std::cout<<" old2_0.5_diff : "<< v_Res.X()<<", "<<v_Res.Y()<<", "<<v_Res.Z()<<std::endl;
   v_Rot = old3_FARO2 - v_Para;
   v_Rot.RotateX(par[3]*TMath::DegToRad());
   v_Rot.RotateY(par[4]*TMath::DegToRad());
   v_Rot.RotateZ(par[5]*TMath::DegToRad());
   v_Res = v_Rot - old3_CAD;
   std::cout<<" old3_0.5_diff : "<< v_Res.X()<<", "<<v_Res.Y()<<", "<<v_Res.Z()<<std::endl;
   v_Rot = old4_FARO2 - v_Para;
   v_Rot.RotateX(par[3]*TMath::DegToRad());
   v_Rot.RotateY(par[4]*TMath::DegToRad());
   v_Rot.RotateZ(par[5]*TMath::DegToRad());
   v_Res = v_Rot - old4_CAD;
   std::cout<<" old4_0.5_diff : "<< v_Res.X()<<", "<<v_Res.Y()<<", "<<v_Res.Z()<<std::endl;
}






