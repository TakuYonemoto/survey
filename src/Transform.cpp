#include <iostream>
#include <fstream>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TMinuit.h>
//#include <setting.h>
 /************************** calculate chi^2 **********************                                        npar            -> number of free parameters involved in minimization
gim             -> computed gradient values (optional)
f               -> function value
par             -> array of parameters
flag            -> to switch between several actions of FCN
*******************************************************************/

using namespace std;
const Int_t nevent = 8;
const Int_t nmaxevent = 8;
Int_t npar = 6;
const char *parname[] = {"x-coo","y-coo","z-coo","alpha","beta","gamma"};
Double_t vstart[] = {100,100,100,100,100,100};
//Double_t vstart[] = {-890,2.15,-579,-0.28,-0.079,-0.063};
Double_t par[nevent];
Double_t err[nevent];

//making the vector
TVector3 v_3d;//3D-scanner-coordinates
TVector3 v_trans;//after translation
TVector3 v_rx;//after rotation around x-axis
TVector3 v_ry;//after rotation around y-axis
TVector3 v_rz;//after rotation around z-axis
TVector3 v_cobra;//cobra-3d-coordinate, not cobra-cylindrical-coordinate

vector<Double_t> xvec3D, yvec3D, zvec3D, xvecCOBRA, yvecCOBRA, zvecCOBRA;
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t flag);
void check(Double_t *par);

void Transform()
{
   
//data setup
   ifstream data("data.dat");
   Double_t temp1, temp2, temp3, temp4, temp5, temp6;
   while (data >> temp1 >> temp2 >> temp3 >> temp4 >> temp5 >> temp6){
       xvec3D.push_back(temp1);
       yvec3D.push_back(temp2);
       zvec3D.push_back(temp3);
       xvecCOBRA.push_back(temp4);
       yvecCOBRA.push_back(temp5);
       zvecCOBRA.push_back(temp6);
   }

   for (Int_t i = 0; i < nevent; i++){
      cout << xvec3D.at(i) << "," << xvecCOBRA.at(i) << endl;
   }
//Minuit setup
   TMinuit *mynuit = new TMinuit(npar);
   mynuit -> SetFCN(fcn);
   Double_t arglist[10];
   arglist[0] = 1;
   Int_t ierflg = 1;
   mynuit -> mnexcm("SET ERR", arglist, 1, ierflg);
   
   for (Int_t i = 0; i < npar; i++){
      mynuit -> mnparm(i, parname[i], vstart[i], 10, 0, 0, ierflg);
   }

   arglist[0] = 10000;
   arglist[1] = 0.1;

   mynuit -> mnexcm("MINOS", arglist, 2, ierflg);
   Double_t fedm, errdef, fmin;
   Int_t npari, nparx, istat;
   mynuit -> mnstat(fmin, fedm, errdef, npari, nparx, istat);
   if (nevent != 2){
      Double_t redchi = fmin/(nevent * 3 - npar);
      cout << "redchi=" << redchi << endl;
   }
   Double_t par_result[nevent];
   for (Int_t i = 0; i < npar; i++){
      Double_t current, error;
      mynuit -> GetParameter(i, current, error);
      par_result[i] = current;
      cout << i << " " << current << " " << error << endl;
   }
   cout << "=========================================================" <<endl;
   check(par_result);
}

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t flag)
{
   Double_t sum = 0;
   Double_t dis2[nevent];

   for (Int_t i = 0; i < nevent; i++){
//      cout << par[0] << "," << par[1] << "," << par[3] << "," << par[4] << "," << par[5] << endl;
      v_3d = TVector3(xvec3D.at(i), yvec3D.at(i), zvec3D.at(i));
      TVector3 v_trans(par[0], par[1], par[2]);
      v_rz = v_3d;
      v_rz.RotateZ(par[5]*TMath::DegToRad());
      v_ry = v_rz;
      v_ry.RotateY(par[4]*TMath::DegToRad());
      v_rx = v_ry;
      v_rx.RotateX(par[3]*TMath::DegToRad());
      v_cobra = v_rx + v_trans;

      Double_t dis_temp1, dis_temp2, dis_temp3;
      dis_temp1 = (v_cobra.X() - xvecCOBRA.at(i));
      dis_temp2 = (v_cobra.Y() - yvecCOBRA.at(i));
      dis_temp3 = (v_cobra.Z() - zvecCOBRA.at(i));
      dis2[i] = dis_temp1 * dis_temp1 + dis_temp2 * dis_temp2 + dis_temp3 * dis_temp3;

      sum += dis2[i];
   }
   f = sum;

}

void check(Double_t *par)
{

   Double_t sum = 0;
   Double_t dis2[nmaxevent];
   for (Int_t i = 0; i < nmaxevent; i++){
//      cout << par[0] << "," << par[1] << "," << par[3] << "," << par[4] << "," << par[5] << endl;
      v_3d = TVector3(xvec3D.at(i), yvec3D.at(i), zvec3D.at(i));
      TVector3 v_trans(par[0], par[1], par[2]);
      v_rz = v_3d;
      v_rz.RotateZ(par[5]*TMath::DegToRad());
      v_ry = v_rz;
      v_ry.RotateY(par[4]*TMath::DegToRad());
      v_rx = v_ry;
      v_rx.RotateX(par[3]*TMath::DegToRad());
      v_cobra = v_rx + v_trans;

      Double_t dis_temp1, dis_temp2, dis_temp3;
      dis_temp1 = (v_cobra.X() - xvecCOBRA.at(i));
      dis_temp2 = (v_cobra.Y() - yvecCOBRA.at(i));
      dis_temp3 = (v_cobra.Z() - zvecCOBRA.at(i));
      dis2[i] = dis_temp1 * dis_temp1 + dis_temp2 * dis_temp2 + dis_temp3 * dis_temp3;
      
      cout << i << ":" << dis_temp1 << "," << dis_temp2 << "," << dis_temp3 <<endl;
      sum += TMath::Sqrt(dis2[i]);
   }
   cout << sum << endl;
}
