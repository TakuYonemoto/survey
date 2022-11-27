#include <iostream>
#include <fstream>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TMath.h>
#include <TMinuit.h>

using namespace TMath;

void LossFunc(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t flag);
void LossFunc_CC(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t flag);
void check(Double_t *par);
void check2(Double_t *par, Double_t *par2);

//Markers
const TVector3 old1_CAD = TVector3(-385.94, -51.184, 1218.4);   //Marker A
const TVector3 old2_CAD = TVector3(-354.39, -45.206, 1092);     //Marker C

const Int_t npoint = 10;
TVector3 v_FARO[npoint], v_COBRA[npoint];

const TVector3 old1_FARO = TVector3(-750.9342, 26.1709, 59.0732);
const TVector3 old2_FARO = TVector3(-734.0348, 26.2725, -51.1029);

const TVector3 old1_COBRA = TVector3(-386.495, -50.981, 1203.715);
const TVector3 old2_COBRA = TVector3(-368.182, -47.530, 1095.251);

Int_t npara = 6;

void FARO_COBRA_calib2018_DS()
{
   TVector3 v_pre = old1_COBRA - old1_CAD;
   std::cout<<" old1_1.5_pre : " << v_pre.X() << ", " << v_pre.Y() << ", " << v_pre.Z() <<std::endl;
   v_pre = old2_COBRA - old2_CAD;
   std::cout<<" old2_1.5_pre : " << v_pre.X() << ", " << v_pre.Y() << ", " << v_pre.Z() <<std::endl;

   v_FARO[0] = TVector3(-728.4323, 0.0554, 0.0046);
   v_FARO[1] = TVector3(-730.7838, -0.3373, -124.4414);
   v_FARO[2] = TVector3(-740.1744, -0.6011, -720.0072);
   v_FARO[3] = TVector3(5.5466, 0.062, -0.008);
   v_FARO[4] = TVector3(4.006, -0.1212, -141.904);
   v_FARO[5] = TVector3(-728.9222, -53.6067, 20.8844);
   v_FARO[6] = TVector3(5.1348, -57.0418, 21.9815);
   v_FARO[7] = TVector3(-524.1375, -326.6278, 146.394);
   v_FARO[8] = TVector3(-395.4759, -363.2516, 144.9649);
   v_FARO[9] = TVector3(-143.3532, -293.5763, 140.5562);
   v_COBRA[0] = TVector3(-359.042, -72.610, 1144.958);
   v_COBRA[1] = TVector3(-359.374, -72.871, 1020.526);
   v_COBRA[2] = TVector3(-359.253, -72.057, 424.711);
   v_COBRA[3] = TVector3(364.576, 53.537, 1155.122);
   v_COBRA[4] = TVector3(365.249, 53.322, 1013.204);
   v_COBRA[5] = TVector3(-350.647, -125.715, 1165.704);
   v_COBRA[6] = TVector3(373.464, -3.133, 1177.149);
   v_COBRA[7] = TVector3(-104.151, -359.632, 1293.950);
   v_COBRA[8] = TVector3(29.011, -373.527, 1294.448);
   v_COBRA[9] = TVector3(265.420, -261.231, 1293.747);

// Fitting
   TMinuit *minuit = new TMinuit(npara);
   minuit->SetFCN(LossFunc);
   Int_t ierflg = 0;
   const char *parname[] = {"x_parallel","y_parallel","z_parallel",
                            "theta_x","theta_y","theta_z"};
   Double_t vstart[] = {0.,0.,0.,0.,0.,0.,0.};
   
   for(Int_t i=0; i < npara; i++){
      minuit->mnparm(i, parname[i],vstart[i], 0.001, 0, 0, ierflg);
   }

   Double_t arglist[10];
   Double_t one_sigma_CL = 0.682689492; // 1sigma ~ 68%
   arglist[0] = TMath::ChisquareQuantile(one_sigma_CL,npara);
   minuit->mnexcm("SET ERR", arglist, 1, ierflg);

   arglist[0] = 10000;
   arglist[1] = 1;
   minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
   Double_t fedm, errdef, fmin;
   Int_t npari, nparx, istat;
   minuit->mnstat(fmin,fedm,errdef,npari,nparx,istat);

   //Double_t redchi = fmin/(4*3 -6);
   std::cout << "fmin = " << fmin << std::endl;
   Double_t par_result[npara], par_err[npara];
   for(Int_t i=0; i<npara; i++){
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

   TMinuit *minuit2 = new TMinuit(npara);
   minuit2->SetFCN(LossFunc_CC);
   for(Int_t i=0; i < npara; i++){
      minuit2->mnparm(i, parname[i],vstart[i], 0.001, 0, 0, ierflg);
   }
   arglist[0] = TMath::ChisquareQuantile(one_sigma_CL,npara);
   minuit2->mnexcm("SET ERR", arglist, 1, ierflg);
   arglist[0] = 10000;
   minuit2->mnexcm("MIGRAD", arglist, 2, ierflg);
   minuit2->mnstat(fmin,fedm,errdef,npari,nparx,istat);
   std::cout << "fmin_CC = " << fmin << std::endl;
   Double_t par_result2[npara], par_err2[npara];
   for(Int_t i=0; i<npara; i++){
      Double_t current, error;
      minuit2->GetParameter(i, current, error);
      par_result2[i] = current;
      par_err2[i] = error;
      if(i==0) std::cout << "parname : value[mm] : error[mm] " << std::endl;
      else if(i==3) std::cout << "parname : value[degree] : error[degree] " << std::endl;
      std::cout << parname[i] << " : " << current << " : " << error << std::endl;
   }
   check2(par_result, par_result2);


   std::ofstream output_file("outputfiles/FARO_COBRA_calib2018_DS.csv");
   for(Int_t i=0; i < npara; i++){
      output_file << par_result[i] << "," << par_err[i] << std::endl;
   }
   output_file.close();
}


void LossFunc(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t flag)
{
   Double_t L = 0.;

   TVector3 v_Rot,v_Res;
   TVector3 v_Para(par[0], par[1], par[2]);
   
   for(Int_t ipoint = 0; ipoint < npoint; ipoint++){   
      v_Rot = v_FARO[ipoint] - v_Para;
      v_Rot.RotateX(par[3]*TMath::DegToRad());
      v_Rot.RotateY(par[4]*TMath::DegToRad());
      v_Rot.RotateZ(par[5]*TMath::DegToRad());
      v_Res = v_COBRA[ipoint] - v_Rot;
      L += v_Res.Mag2();
   }

   v_Rot = old1_FARO - v_Para;
   v_Rot.RotateX(par[3]*TMath::DegToRad());
   v_Rot.RotateY(par[4]*TMath::DegToRad());
   v_Rot.RotateZ(par[5]*TMath::DegToRad());
   v_Res = old1_COBRA - v_Rot;
   L += v_Res.Mag2();

   v_Rot = old2_FARO - v_Para;
   v_Rot.RotateX(par[3]*TMath::DegToRad());
   v_Rot.RotateY(par[4]*TMath::DegToRad());
   v_Rot.RotateZ(par[5]*TMath::DegToRad());
   v_Res = old2_COBRA - v_Rot;
   L += v_Res.Mag2();

   f = L;

}

void LossFunc_CC(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t flag)
{

   Double_t L = 0.;

   TVector3 v_Rot,v_Res;
   TVector3 v_Para(par[0], par[1], par[2]);

   v_Rot = old1_COBRA - v_Para;
   v_Rot.RotateX(par[3]*TMath::DegToRad());
   v_Rot.RotateY(par[4]*TMath::DegToRad());
   v_Rot.RotateZ(par[5]*TMath::DegToRad());
   v_Res = old1_CAD - v_Rot;
   L += v_Res.Mag2();
   
   v_Rot = old2_COBRA - v_Para;
   v_Rot.RotateX(par[3]*TMath::DegToRad());
   v_Rot.RotateY(par[4]*TMath::DegToRad());
   v_Rot.RotateZ(par[5]*TMath::DegToRad());
   v_Res = old2_CAD - v_Rot;
   L += v_Res.Mag2();

   f = L;

}

void check(Double_t *par)
{
   Double_t L = 0;

   TVector3 v_Rot,v_Res;
   TVector3 v_Para(par[0], par[1], par[2]);
   for(Int_t ipoint=0; ipoint < npoint; ipoint++){
      v_Rot = v_FARO[ipoint] - v_Para;
      v_Rot.RotateX(par[3]*TMath::DegToRad());
      v_Rot.RotateY(par[4]*TMath::DegToRad());
      v_Rot.RotateZ(par[5]*TMath::DegToRad());
      v_Res = v_COBRA[ipoint] - v_Rot;
      std::cout << ipoint << " : " << v_Res.X() << ", " << v_Res.Y() << ", " << v_Res.Z() << std::endl;
      L += v_Res.Mag2();
   }
 
   std::cout << "LossFunc = " << L << endl;
   std::cout<<"    ID    : x_diff, y_diff, z_diff ([mm] in COBRA coodinate system)"<<std::endl;
   
   v_Rot = old1_FARO - v_Para;
   v_Rot.RotateX(par[3]*TMath::DegToRad());
   v_Rot.RotateY(par[4]*TMath::DegToRad());
   v_Rot.RotateZ(par[5]*TMath::DegToRad());
   v_Res = old1_COBRA - v_Rot;
   L += v_Res.Mag2();
   std::cout<<" old1_1.5 : "<< v_Res.X()<<", "<<v_Res.Y()<<", "<<v_Res.Z()<<std::endl;
   
   v_Rot = old2_FARO - v_Para;
   v_Rot.RotateX(par[3]*TMath::DegToRad());
   v_Rot.RotateY(par[4]*TMath::DegToRad());
   v_Rot.RotateZ(par[5]*TMath::DegToRad());
   v_Res = old2_COBRA - v_Rot;
   L += v_Res.Mag2();
   std::cout<<" old2_1.5 : "<< v_Res.X()<<", "<<v_Res.Y()<<", "<<v_Res.Z()<<std::endl;

   std::cout << "LossFunc_1.5 = " << L << endl;
   std::cout<<"***************************************************"<<std::endl;
}

void check2(Double_t *par, Double_t *par2)
{
   Double_t L = 0;

   TVector3 v_Rot,v_Res;
   TVector3 v_Para(par[0], par[1], par[2]);
   TVector3 v_Para2(par2[0], par2[1], par2[2]);
 
   std::cout<<"    ID    : x_diff, y_diff, z_diff ([mm] in COBRA coodinate system)"<<std::endl;
   
   v_Rot = old1_FARO - v_Para;
   v_Rot.RotateX(par[3]*TMath::DegToRad());
   v_Rot.RotateY(par[4]*TMath::DegToRad());
   v_Rot.RotateZ(par[5]*TMath::DegToRad());
   v_Res = old1_COBRA - v_Rot;
   std::cout<<" old1_1.5 : "<< v_Res.X()<<", "<<v_Res.Y()<<", "<<v_Res.Z()<<std::endl;
   v_Rot = v_Rot - v_Para2;
   v_Rot.RotateX(par2[3]*TMath::DegToRad());
   v_Rot.RotateY(par2[4]*TMath::DegToRad());
   v_Rot.RotateZ(par2[5]*TMath::DegToRad());
   v_Res = old1_COBRA - v_Rot;
   L += v_Res.Mag2();
   std::cout<<" old1_1.5 : "<< v_Res.X()<<", "<<v_Res.Y()<<", "<<v_Res.Z()<<std::endl;
   
   v_Rot = old2_FARO - v_Para;
   v_Rot.RotateX(par[3]*TMath::DegToRad());
   v_Rot.RotateY(par[4]*TMath::DegToRad());
   v_Rot.RotateZ(par[5]*TMath::DegToRad());
   v_Res = old1_COBRA - v_Rot;
   std::cout<<" old1_1.5 : "<< v_Res.X()<<", "<<v_Res.Y()<<", "<<v_Res.Z()<<std::endl;
   v_Rot = v_Rot - v_Para2;
   v_Rot.RotateX(par2[3]*TMath::DegToRad());
   v_Rot.RotateY(par2[4]*TMath::DegToRad());
   v_Rot.RotateZ(par2[5]*TMath::DegToRad());
   v_Res = old2_COBRA - v_Rot;
   L += v_Res.Mag2();
   std::cout<<" old2_1.5 : "<< v_Res.X()<<", "<<v_Res.Y()<<", "<<v_Res.Z()<<std::endl;

   std::cout << "LossFunc2 = " << L << endl;
   std::cout<<"***************************************************"<<std::endl;

}

