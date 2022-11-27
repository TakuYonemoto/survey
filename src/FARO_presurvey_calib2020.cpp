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

std::vector<TVector3> V,V_pre; 
std::vector<Double_t> ID,ID_pre; 
Int_t skip_rows = 2;
Int_t ID_size,ID_pre_size;
Int_t npar = 6;

void FARO_presurvey_calib2020(){
   
   Double_t x,y,z;
   Int_t id;
   TVector3 v;
   std::string str;
   
   std::ifstream ifs("data/references_survey_2019.csv");
   if(!ifs){
      cout << "Error: File cannot open!" << endl;
   } 
   for(Int_t i=0; i<skip_rows; i++){
      getline(ifs,str);
   }
   while(getline(ifs,str)){
      sscanf(str.data(), "%d,%lf,%lf,%lf",&id,&x,&y,&z);
      v = TVector3(x,y,z);
      V.push_back(v);
      ID.push_back(id);
   }
   ifs.close();

   ifs.open("data/references_presurvey_2020.csv");
   if(!ifs){
      cout << "Error: File cannot open!" << endl;
   } 
   for(Int_t i=0; i<skip_rows; i++){
      getline(ifs,str);
   }
   while(getline(ifs,str)){
      sscanf(str.data(), "%d,%lf,%lf,%lf",&id,&x,&y,&z);
      v = TVector3(x,y,z);
      V_pre.push_back(v);
      ID_pre.push_back(id);
   }
   ifs.close();

   ID_size = ID.size();
   ID_pre_size = ID_pre.size();

   
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
   std::ofstream output_file("outputfiles/FARO_presurvey_calib2020.csv");
   for(Int_t i=0; i < npar; i++){
      output_file << par_result[i] << "," << par_err[i] << std::endl;
   }
   output_file.close();

   TVector3 v_Rot;
   TVector3 v_Para(par_result[0], par_result[1], par_result[2]);
   output_file.open("outputfiles/FARO_survey_virtual_2020.csv");
   output_file << "reference_id,x,y,z" << std::endl;
   output_file << "(),(mm),(mm),(mm)" << std::endl;
   for(Int_t i=0; i < ID_size; i++){
      output_file << ID[i] << "," << V[i].X() << "," << V[i].Y() << "," << V[i].Z() << std::endl;
   }
   for(Int_t i=0; i < ID_pre_size; i++){
      if(ID_pre[i] < -200){
         v_Rot = V_pre[i] - v_Para;
         v_Rot.RotateX(par_result[3]*TMath::DegToRad());
         v_Rot.RotateY(par_result[4]*TMath::DegToRad());
         v_Rot.RotateZ(par_result[5]*TMath::DegToRad());
         output_file << ID_pre[i] << "," << v_Rot.X() << "," << v_Rot.Y() << "," << v_Rot.Z() << std::endl;
      }else{
         continue;
      }
   }
   output_file.close();
}


void LossFunc(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t flag)
{
   Double_t L = 0.;

   Int_t id_1, id_2, j_match;
   TVector3 v_Rot,v_Res;
   TVector3 v_Para(par[0], par[1], par[2]);
   for(Int_t i=0; i < ID_size; i++){
      id_1 = ID[i];
      for(Int_t j=0; j < ID_pre_size; j++){
         id_2 = ID_pre[j];
         if(id_2 == id_1){
            j_match = j;
         }
      }
      v_Rot = V_pre[j_match] - v_Para;
      v_Rot.RotateX(par[3]*TMath::DegToRad());
      v_Rot.RotateY(par[4]*TMath::DegToRad());
      v_Rot.RotateZ(par[5]*TMath::DegToRad());
      v_Res = V[i] - v_Rot;
      L += v_Res.Mag2();
   }
   f = L;

}

void check(Double_t *par)
{
   Double_t L = 0;
   
   Int_t id_1, id_2, j_match;
   TVector3 v_Rot,v_Res;
   TVector3 v_Para(par[0], par[1], par[2]);
   for(Int_t i=0; i < ID_size; i++){
      id_1 = ID[i];
      for(Int_t j=0; j < ID_pre_size; j++){
         id_2 = ID_pre[j];
         if(id_2 == id_1){
            j_match = j;
         }
      }
      v_Rot = V_pre[j_match] - v_Para;
      v_Rot.RotateX(par[3]*TMath::DegToRad());
      v_Rot.RotateY(par[4]*TMath::DegToRad());
      v_Rot.RotateZ(par[5]*TMath::DegToRad());
      v_Res = V[i] - v_Rot;
      std::cout << " " << id_1 << " : " << v_Res.X() << ", " << v_Res.Y() << ", " << v_Res.Z() << std::endl;
      L += v_Res.Mag2();
   }

   std::cout << "LossFunc = " << L << endl;
 
}






