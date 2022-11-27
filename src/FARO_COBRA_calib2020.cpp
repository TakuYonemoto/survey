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

//Markers
TVector3 old1_CAD = TVector3( 380.17,  83.901, -1218.4);   //Marker A
TVector3 old3_CAD = TVector3(-385.94, -51.184, -1218.4); //Marker B
TVector3 old2_CAD = TVector3( 348.51,  78.733, -1092);     //Marker C
TVector3 old4_CAD = TVector3(-354.42, -45.213, -1092);   //Marker D

std::vector<TVector3> V_FARO,V_COBRA,Verr;
std::vector<Double_t> ID_FARO,ID_COBRA; 
Int_t skip_rows = 2;
Int_t ID_FARO_size,ID_COBRA_size;
Int_t npar = 6;
Int_t nevent = 4;

void FARO_COBRA_calib2020(){
   
   Double_t x,y,z,dx,dy,dz;
   Int_t id;
   TVector3 v;
   std::string str;
   
   std::ifstream ifs("outputfiles/references_survey_virtual_2020.csv");
   if(!ifs){
      cout << "Error: File ifs(FARO) cannot open!" << endl;
   } 
   for(Int_t i=0; i<skip_rows; i++){
      getline(ifs,str);
   }
   while(getline(ifs,str)){
      sscanf(str.data(), "%d,%lf,%lf,%lf",&id,&x,&y,&z);
      v = TVector3(x,y,z);
      V_FARO.push_back(v);
      ID_FARO.push_back(id);
   }
   ifs.close();

   ifs.open("data/references_cobra_2020.csv");
   if(!ifs){
      cout << "Error: File ifs(COBRA) cannot open!" << endl;
   } 
   for(Int_t i=0; i<skip_rows; i++){
      getline(ifs,str);
   }
   while(getline(ifs,str)){
      sscanf(str.data(), "%d,%lf,%lf,%lf,%lf,%lf,%lf",&id,&x,&y,&z,&dx,&dy,&dz);
      v = TVector3(x,y,z);
      V_COBRA.push_back(v);
      ID_COBRA.push_back(id);
      v = TVector3(dx,dy,dz);
      Verr.push_back(v);
   }
   ifs.close();

   ID_FARO_size = ID_FARO.size();
   ID_COBRA_size = ID_COBRA.size();

   TMinuit *minuit = new TMinuit(npar);
   minuit->SetFCN(LossFunc);
   Int_t ierflg = 0;
   const char *parname[] = {"x_parallel","y_parallel","z_parallel","theta_x","theta_y","theta_z"};
   Double_t vstart[] = {375.,-9.,-1145.,359,180,190};
   
   for(Int_t i=0; i < npar; i++){
      minuit->mnparm(i, parname[i],vstart[i], 0.1, 0, 0, ierflg);
   }

   Double_t arglist[10];
   arglist[0] = TMath::ChisquareQuantile(0.682689492,1);
   //arglist[0] = TMath::ChisquareQuantile(0.682689492,npar);
   minuit->mnexcm("SET ERR", arglist, 1, ierflg);

   arglist[0] = 10000;
   arglist[1] = 1;
   minuit->mnexcm("MINOS", arglist, 2, ierflg);
   Double_t fedm, errdef, fmin;
   Int_t npari, nparx, istat;
   minuit->mnstat(fmin,fedm,errdef,npari,nparx,istat);

   //Double_t redchi = fmin/(nevent*3 - 6);
   std::cout << "LossFunc = " << fmin << std::endl;
   if(nevent != 2) std::cout << "redchi = " << fmin/(nevent - npar) << std::endl;
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
   std::ofstream output_file("outputfiles/FARO_COBRA_calib2020.csv");
   for(Int_t i=0; i < npar; i++){
      output_file << par_result[i] << "," << par_err[i] << std::endl;
   }
   output_file.close();


   TVector3 v_Rot;
   TVector3 v_Para(par_result[0], par_result[1], par_result[2]);
   output_file.open("outputfiles/references_cobra_virtual_2020.csv");
   output_file << "reference_id,x,y,z" << std::endl;
   output_file << "(),(mm),(mm),(mm)" << std::endl;
   for(Int_t i=0; i < ID_FARO_size; i++){
      v_Rot = V_FARO[i] - v_Para;
      v_Rot.RotateX(par_result[3]*TMath::DegToRad());
      v_Rot.RotateY(par_result[4]*TMath::DegToRad());
      v_Rot.RotateZ(par_result[5]*TMath::DegToRad());
      output_file << ID_FARO[i] << "," << v_Rot.X() << "," << v_Rot.Y() << "," << v_Rot.Z() << std::endl;
   }
   output_file.close();
}


void LossFunc(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t flag)
{
   Double_t L = 0.;

   TVector3 v_Rot,v_Res;
   TVector3 v_Para(par[0], par[1], par[2]);
   for(Int_t i=0; i < ID_FARO_size; i++){
      for(Int_t j=0; j < ID_COBRA_size; j++){
         //if(ID_COBRA[j] == ID_FARO[i]){
         if(ID_COBRA[j] == ID_FARO[i] && ID_COBRA[j] > -200){
            v_Rot = V_FARO[i] - v_Para;
            v_Rot.RotateX(par[3]*TMath::DegToRad());
            v_Rot.RotateY(par[4]*TMath::DegToRad());
            v_Rot.RotateZ(par[5]*TMath::DegToRad());
            v_Res = V_COBRA[j] - v_Rot;
            //L += v_Res.Mag2();
            L += v_Res.Mag2()/Verr[j].Mag2();
         }else{
            continue;
         }
      }
   }
   f = L;

}

void check(Double_t *par)
{
   Double_t L = 0;
   
   TVector3 v_Rot,v_Res;
   TVector3 v_Para(par[0], par[1], par[2]);
   for(Int_t i=0; i < ID_FARO_size; i++){
      for(Int_t j=0; j < ID_COBRA_size; j++){
         //if(ID_COBRA[j] == ID_FARO[i]){
         if(ID_COBRA[j] == ID_FARO[i] && ID_COBRA[j] > -200){
            v_Rot = V_FARO[i] - v_Para;
            v_Rot.RotateX(par[3]*TMath::DegToRad());
            v_Rot.RotateY(par[4]*TMath::DegToRad());
            v_Rot.RotateZ(par[5]*TMath::DegToRad());
            v_Res = V_COBRA[j] - v_Rot;
            std::cout << ID_FARO[i] << " : " << v_Res.X() << ", " << v_Res.Y() << ", " << v_Res.Z() << std::endl;
            //std::cout << ID_1[i] << "_err:" << Verr[j].X() << ", " << Verr[j].Y() << ", " << Verr[j].Z() << std::endl;
            //L += v_Res.Mag2();
            L += v_Res.Mag2()/Verr[j].Mag2();
         }else{
            continue;
         }
      }
   }

   std::cout << "LossFunc = " << L << endl;
   std::cout << "redchi = " << L/(nevent - npar) << endl;
   //std::cout << "redchi = " << L/(nevent*3 - npar) << endl;
 
}

