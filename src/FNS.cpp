#include "TROOT.h"
#include "TBrowser.h"
#include <sstream>
#include <stdio.h>
#include <vector>
#include <fstream>
#include <TVector3.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TMatrix.h>

using namespace TMath;
void Dis(TVector3 a, TVector3 &a_prime, Double_t &r);
void CalcML(TVector3 a[], TVectorD q, TMatrixD R, TMatrixD W, TMatrixD &M, TMatrixD &L)

const Double_t xmax_SidePlane = 18.1544845 * 2;
const Double_t xmin_SidePlane = -0.84551545;    //x2
const Double_t x_TopPlane     = -20.;           //x1
const Double_t y_SidePlane    = 5.4923668;      //y2
const Double_t y_TopPlane     = 3.1;            //y1
const Double_t z_PCB          = 61.65;
const Double_t slope_approx   = -0.12489852;    // (y2-y1)/(x2-x1) : approximate side curve by line
const Double_t const_approx   = 7.9903372;      // slope_approx * x1 + y1 : const of line
const Double_t slope_square   = 0.015599641;    // slope_approx^2

const Int_t NRepeat = 10;

void FNS(TVector3 a[], TMatrixD &R){

   TMatrixD T[3];
   for(Int_t i=0; i<3; i++) {T[i].ResizeTo(4,6);}
   T[0](0,0) = T[0](2,2) = T[0](2,5) = T[1](0,1) = T[1](3,0) = T[1](3,3) = T[2](0,2) = T[2](1,1) = T[2](1,4) = -1.;
   T[0](0,3) = T[0](3,1) = T[0](3,4) = T[1](0,4) = T[1](1,2) = T[1](1,5) = T[2](0,5) = T[2](2,0) = T[2](2,3) = 1.;

   TMatrixD V0(6,6);
   V0(0,0) = V0(1,1) = V0(2,2) = V0(3,3) = V0(4,4) = V0(5,5) = 1.; //共分散行列を単位行列で代用(≒異方性誤差無しを仮定)

   TMatrixD V0_xi[3][3];
   for(Int_t k=0; k<3; k++){
      for(Int_t l=0; l<3; l++){
         V0_xi[k][l].ResizeTo(4,4);
         V0_xi[k][l] =  T[k] * V0 * T[l];
      }
   }

   TVectorD q,q0; q.ResizeTo(4); q0.ResizeTo(4);
   q = (0,0,0,0); q0 = (0,0,0,0);
   TMatrixD R(3,3); R(0,0) = R(1,1) = R(2,2) = 1.;

   TMatrixD X,M,L,W;
   W.ResizeTo(3,3); W(0,0) = W(1,1) = W(2,2) = 1.;
   M.ResizeTo(4,4); L.ResizeTo(4,4); X.ResizeTo(4,4);

   CalcML(a,q,R,W,M,L);

   Double_t lambda,lambda_1;

   for(Int_t i=0; i<NRepeat; i++){

      X = M - L;
      auto XEigen = new TMatrixDEigen(X);
      lambda = XEigen.GetEigenValues()(3,3);
      lambda_1 = XEigen.GetEigenValues()(2,2);
      if(lambda < lambda_1){
         q = TMatrixDColumn(XEigen.GetEigenVectors(),3);
         q = q * (1./sqrt(q.Norm2Sqr()));
      }else{
         q = TMatrixDColumn(XEigen.GetEigenVectors(),2);
         q = q * (1./sqrt(q.Norm2Sqr()));
      }

      if(sqrt((q-q0).Norm2Sqr < 0.1) || sqrt((q+q0).Norm2Sqr) < 0.1){
         continue;
      }else if{
         R(0,0) = q(0)*q(0)+q(1)*q(1)-q(2)*q(2)-q(3)*q(3);
         R(0,1) = (q(1)*q(2)-q(0)*q(3))*2;
         R(0,2) = (q(1)*q(3)+q(0)*q(2))*2;
         R(1,0) = (q(1)*q(2)+q(0)*q(3))*2;
         R(1,1) = q(0)*q(0)-q(1)*q(1)+q(2)*q(2)-q(3)*q(3);
         R(1,2) = (q(2)*q(3)-q(0)*q(1))*2;
         R(2,0) = (q(1)*q(3)-q(0)*q(2))*2;
         R(2,1) = (q(2)*q(3)+q(0)*q(1))*2;
         R(2,2) = q(0)*q(0)-q(1)*q(1)-q(2)*q(2)+q(3)*q(3);

         q0 = q;

         for(Int_t k=0; k<3; k++){
            for(Int_t l=0; l<3; l++){
               W(k,l) = q * (V0_xi[k][l] * q);
            }
         }
         W.Invert();
         CalcML(a,q,R,W,M,L);
      }
   }
   std::cout << "lambda = " << lambda << std::endl;
}

void Dis(TVector3 a, TVector3 &a_prime, Double_t &r){
   Double_t x = a.X();
   Double_t y = a.Y();
   Double_t z = a.Z();
   Double_t m = 0;
   Int_t min_num = 0;
   Double_t R[9];

   Double_t xS1 = Power(x - xmax_SidePlane, 2);
   Double_t xS2 = Power(x - xmin_SidePlane, 2);
   Double_t xT = Power(x - x_TopPlane, 2);

   Double_t yS1 = Power(y - y_SidePlane,2);
   Double_t yS2 = Power(y + y_SidePlane,2);
   Double_t yT = Power(Abs(y) - y_TopPlane, 2);

   Double_t zS1 = Power(z - z_PCB,2);
   Double_t zS2 = Power(z + z_PCB,2);

   Double_t x_th1 = (x + slope_approx * (y - const_approx))/(slope_square + 1);
   Double_t x_th2 = (x - slope_approx * (y + const_approx))/(slope_square + 1);
   Double_t y_th1 = slope_approx * x_th1 + const_approx;
   Double_t y_th2 = -slope_approx * x_th2 - const_approx;
   Double_t d_approx1 = Power((slope_approx * x - y + const_approx),2)/(slope_square + 1);
   Double_t d_approx2 = Power((-slope_approx * x - y - const_approx),2)/(slope_square + 1);
   Double_t y_th = Abs(slope_approx * x + const_approx);

   //平面・大1(y>0)
   if     (x >= xmax_SidePlane && Abs(z) >= z_PCB)       {R[0] = Sqrt(xS1 + yS1 + Min(zS1,zS2));}
   else if(x <= xmin_SidePlane && Abs(z) >= z_PCB)       {R[0] = Sqrt(xS2 + yS1 + Min(zS1,zS2));}
   else if(x >= xmax_SidePlane)                          {R[0] = Sqrt(xS1 + yS1);}
   else if(x <= xmin_SidePlane)                          {R[0] = Sqrt(xS2 + yS1);}
   else if(Abs(z) >= z_PCB)                              {R[0] = Sqrt(yS1 + Min(zS1,zS2));}
   else                                                  {R[0] = Sqrt(yS1);}
   //平面・大2(y<0)
   if     (x >= xmax_SidePlane && Abs(z) >= z_PCB)       {R[1] = Sqrt(xS1 + yS2 + Min(zS1,zS2));}
   else if(x <= xmin_SidePlane && Abs(z) >= z_PCB)       {R[1] = Sqrt(xS2 + yS2 + Min(zS1,zS2));}
   else if(x >= xmax_SidePlane)                          {R[1] = Sqrt(xS1 + yS2);}
   else if(x <= xmin_SidePlane)                          {R[1] = Sqrt(xS2 + yS2);}
   else if(Abs(z) >= z_PCB)                              {R[1] = Sqrt(yS2 + Min(zS1,zS2));}
   else                                                  {R[1] = Sqrt(yS2);}
   //平面・小1(z>0)
   if     (x >= xmax_SidePlane && Abs(y) >= y_SidePlane) {R[2] = Sqrt(xS1 + Min(yS1,yS2) + zS1);}
   else if(x <= xmin_S1idePlane && Abs(y) >= y_SidePlane){R[2] = Sqrt(xS2 + Min(yS1,yS2) + zS1);}
   else if(x >= xmax_SidePlane)                          {R[2] = Sqrt(xS1 + zS1);}
   else if(x <= xmin_SidePlane)                          {R[2] = Sqrt(xS2 + zS1);}
   else if(Abs(y) >= y_SidePlane)                        {R[2] = Sqrt(Min(yS1,yS2) + zS1);}
   else                                                  {R[2] = Sqrt(zS1);}
   //平面・小2(z<0)
   if     (x >= xmax_SidePlane && Abs(y) >= y_SidePlane) {R[3] = Sqrt(xS1 + Min(yS1,yS2) + zS2);}
   else if(x <= xmin_SidePlane && Abs(y) >= y_SidePlane) {R[3] = Sqrt(xS2 + Min(yS1,yS2) + zS2);}
   else if(x >= xmax_SidePlane)                          {R[3] = Sqrt(xS1 + zS2);}
   else if(x <= xmin_SidePlane)                          {R[3] = Sqrt(xS2 + zS2);}
   else if(Abs(y) >= y_SidePlane)                        {R[3] = Sqrt(Min(yS1,yS2) + zS2);}
   else                                                  {R[3] = Sqrt(zS2);}
   //上面
   if     (Abs(y) > y_TopPlane && Abs(z) > z_PCB)        {R[4] = Sqrt(xT + yT + Min(zS1,zS2));}
   else if(Abs(y) > y_TopPlane)                          {R[4] = Sqrt(xT + yT);}
   else if(Abs(z) > z_PCB)                               {R[4] = Sqrt(xT + Min(zS1,zS2));}
   else                                                  {R[4] = Sqrt(xT);}
   //曲面・大1(y>0)
   if     (x_th1 >= xmin_SidePlane && Abs(z) >= z_PCB)   {R[5] = Sqrt(xS2 + yS1 + Min(zS1,zS2));}
   else if(x_th1 <= x_TopPlane && Abs(z) >= z_PCB)       {R[5] = Sqrt(xT + yT + Min(zS1,zS2));}
   else if(x_th1 >= xmin_SidePlane)                      {R[5] = Sqrt(xS2 + yS1);}
   else if(x_th1 <= x_TopPlane)                          {R[5] = Sqrt(xT + yT);}
   else if(Abs(z) >= z_PCB)                              {R[5] = Sqrt(d_approx1 + Min(zS1,zS2));}
   else                                                  {R[5] = Sqrt(d_approx1);}
   //曲面・大2(y<0)
   if     (x_th2 >= xmin_SidePlane && Abs(z) >= z_PCB)   {R[6] = Sqrt(xS2 + yS2 + Min(zS1,zS2));}
   else if(x_th2 <= x_TopPlane && Abs(z) >= z_PCB)       {R[6] = Sqrt(xT + yT + Min(zS1,zS2));}
   else if(x_th2 >= xmin_SidePlane)                      {R[6] = Sqrt(xS2 + yS2);}
   else if(x_th2 <= x_TopPlane)                          {R[6] = Sqrt(xT + yT);}
   else if(Abs(z) >= z_PCB)                              {R[6] = Sqrt(d_approx2 + Min(zS1,zS2));}
   else                                                  {R[6] = Sqrt(d_approx2);}
   //曲面・小1(z>0)
   if     (x >= xmin_SidePlane && Abs(y) >= y_SidePlane) {R[7] = Sqrt(xS2 + Min(yS1,yS2) + zS1);}
   else if(x <= x_TopPlane && Abs(y) >= y_TopPlane)      {R[7] = Sqrt(xT + yT + zS1);}
   else if(x >= xmin_SidePlane)                          {R[7] = Sqrt(xS2 + zS1);}
   else if(x <= x_TopPlane)                              {R[7] = Sqrt(xT + zS1);}
   else if(Abs(y) <= y_th)                               {R[7] = Sqrt(zS1);}
   else if(y_th1 >= y_SidePlane)                         {R[7] = Sqrt(xS2 + yS1 + zS1);}
   else if(y_th2 <= -y_SidePlane)                        {R[7] = Sqrt(xS2 + yS2 + zS1);}
   else                                                  {R[7] = Sqrt(Min(d_approx1,d_approx2) + zS1);}
   //曲面・小2(z<0)
   if     (x >= xmin_SidePlane && Abs(y) >= y_SidePlane) {R[8] = Sqrt(xS2 + Min(yS1,yS2) + zS2);}
   else if(x <= x_TopPlane && Abs(y) >= y_TopPlane)      {R[8] = Sqrt(xT + yT + zS2);}
   else if(x >= xmin_SidePlane)                          {R[8] = Sqrt(xS2 + zS2);}
   else if(x <= x_TopPlane)                              {R[8] = Sqrt(xT + zS2);}
   else if(Abs(y) <= y_th)                               {R[8] = Sqrt(zS2);}
   else if(y_th1 >= y_SidePlane)                         {R[8] = Sqrt(xS2 + yS1 + zS2);}
   else if(y_th2 <= -y_SidePlane)                        {R[8] = Sqrt(xS2 + yS2 + zS2);}
   else                                                  {R[8] = Sqrt(Min(d_approx1,d_approx2) + zS2);}
   
   m = R[0];
   for(Int_t i=1; i<9; i++){
      if(m > R[i]){
         m = R[i];
         min_num = i;
      }else{continue;}
   }
   if(min_num == 0){
      if     (x >= xmax_SidePlane && z >= z_PCB)       {Double_t a_p[3] = {xmax_SidePlane, y_SidePlane, z_PCB};}
      else if(x >= xmax_SidePlane && z <= -z_PCB)      {Double_t a_p[3] = {xmax_SidePlane, y_SidePlane, -z_PCB};}
      else if(x <= xmin_SidePlane && z >= z_PCB)       {Double_t a_p[3] = {xmin_SidePlane, y_SidePlane, z_PCB};}
      else if(x <= xmin_SidePlane && z <= -z_PCB)      {Double_t a_p[3] = {xmin_SidePlane, y_SidePlane, -z_PCB};}
      else if(x >= xmax_SidePlane)                     {Double_t a_p[3] = {xmax_SidePlane, y_SidePlane, z};}
      else if(x <= xmin_SidePlane)                     {Double_t a_p[3] = {xmin_SidePlane, y_SidePlane, z};}
      else if(z >= z_PCB)                              {Double_t a_p[3] = {x, y_SidePlane, z_PCB};}
      else if(z <= -z_PCB)                             {Double_t a_p[3] = {x, y_SidePlane, -z_PCB};}
      else                                             {Double_t a_p[3] = {x, y_SidePlane, z};}
   }else if(min_num == 1){
      if     (x >= xmax_SidePlane && z >= z_PCB)       {Double_t a_p[3] = {xmax_SidePlane, -y_SidePlane, z_PCB};}
      else if(x >= xmax_SidePlane && z <= -z_PCB)      {Double_t a_p[3] = {xmax_SidePlane, -y_SidePlane, -z_PCB};}
      else if(x <= xmin_SidePlane && z >= z_PCB)       {Double_t a_p[3] = {xmin_SidePlane, -y_SidePlane, z_PCB};}
      else if(x <= xmin_SidePlane && z <= -z_PCB)      {Double_t a_p[3] = {xmin_SidePlane, -y_SidePlane, -z_PCB};}
      else if(x >= xmax_SidePlane)                     {Double_t a_p[3] = {xmax_SidePlane, -y_SidePlane, z};}
      else if(x <= xmin_SidePlane)                     {Double_t a_p[3] = {xmin_SidePlane, -y_SidePlane, z};}
      else if(z >= z_PCB)                              {Double_t a_p[3] = {x, -y_SidePlane, z_PCB};}
      else if(z <= -z_PCB)                             {Double_t a_p[3] = {x, -y_SidePlane, -z_PCB};}
      else                                             {Double_t a_p[3] = {x, -y_SidePlane, z};}
   }else if(min_num == 2){
      if     (x >= xmax_SidePlane && y >= y_SidePlane) {Double_t a_p[3] = {xmax_SidePlane, y_SidePlane, z_PCB};}
      else if(x >= xmax_SidePlane && y <= -y_SidePlane){Double_t a_p[3] = {xmax_SidePlane, -y_SidePlane, z_PCB};}
      else if(x <= xmin_SidePlane && y >= y_SidePlane) {Double_t a_p[3] = {xmin_SidePlane, y_SidePlane, z_PCB};}
      else if(x <= xmin_SidePlane && y <= -y_SidePlane){Double_t a_p[3] = {xmin_SidePlane, -y_SidePlane, z_PCB};}
      else if(x >= xmax_SidePlane)                     {Double_t a_p[3] = {xmax_SidePlane, y, z_PCB};}
      else if(x <= xmin_SidePlane)                     {Double_t a_p[3] = {xmin_SidePlane, y, z_PCB};}
      else if(y >= y_SidePlane)                        {Double_t a_p[3] = {x, y_SidePlane, z_PCB};}
      else if(y <= -y_SidePlane)                       {Double_t a_p[3] = {x, -y_SidePlane, z_PCB};}
      else                                             {Double_t a_p[3] = {x, y, z_PCB};}
   }else if(min_num == 3){
      if     (x >= xmax_SidePlane && y >= y_SidePlane) {Double_t a_p[3] = {xmax_SidePlane, y_SidePlane, -z_PCB};}
      else if(x >= xmax_SidePlane && y <= -y_SidePlane){Double_t a_p[3] = {xmax_SidePlane, -y_SidePlane, -z_PCB};}
      else if(x <= xmin_SidePlane && y >= y_SidePlane) {Double_t a_p[3] = {xmin_SidePlane, y_SidePlane, -z_PCB};}
      else if(x <= xmin_SidePlane && y <= -y_SidePlane){Double_t a_p[3] = {xmin_SidePlane, -y_SidePlane, -z_PCB};}
      else if(x >= xmax_SidePlane)                     {Double_t a_p[3] = {xmax_SidePlane, y, -z_PCB};}
      else if(x <= xmin_SidePlane)                     {Double_t a_p[3] = {xmin_SidePlane, y, -z_PCB};}
      else if(y >= y_SidePlane)                        {Double_t a_p[3] = {x, y_SidePlane, -z_PCB};}
      else if(y <= -y_SidePlane)                       {Double_t a_p[3] = {x, -y_SidePlane, -z_PCB};}
      else                                             {Double_t a_p[3] = {x, y, -z_PCB};}
   }else if(min_num ==4){
      if     (y > y_TopPlane && z > z_PCB)             {Double_t a_p[3] = {x_TopPlane, y_TopPlane, z_PCB};}
      else if(y < -y_TopPlane && z > z_PCB)            {Double_t a_p[3] = {x_TopPlane, -y_TopPlane, z_PCB};}
      else if(y > y_TopPlane && z < -z_PCB)            {Double_t a_p[3] = {x_TopPlane, y_TopPlane, -z_PCB};}
      else if(y < -y_TopPlane && z < -z_PCB)           {Double_t a_p[3] = {x_TopPlane, -y_TopPlane, -z_PCB};}
      else if(y > y_TopPlane)                          {Double_t a_p[3] = {x_TopPlane, y_TopPlane, z};}
      else if(y < -y_TopPlane)                         {Double_t a_p[3] = {x_TopPlane, -y_TopPlane, z};}
      else if(z > z_PCB)                               {Double_t a_p[3] = {x_TopPlane, y, z_PCB};}
      else if(z < -z_PCB)                              {Double_t a_p[3] = {x_TopPlane, y, -z_PCB};}
      else                                             {Double_t a_p[3] = {x_TopPlane, y, z};}
   }else if(min_num == 5){
      if     (x_th1 >= xmin_SidePlane && z >= z_PCB)   {Double_t a_p[3] = {xmin_SidePlane, y_SidePlane, z_PCB};}
      else if(x_th1 >= xmin_SidePlane && z <= -z_PCB)  {Double_t a_p[3] = {xmin_SidePlane, y_SidePlane, -z_PCB};}
      else if(x_th1 <= x_TopPlane && z >= z_PCB)       {Double_t a_p[3] = {x_TopPlane, y_TopPlane, z_PCB};}
      else if(x_th1 <= x_TopPlane && z <= -z_PCB)      {Double_t a_p[3] = {x_TopPlane, y_TopPlane, -z_PCB};}
      else if(x_th1 >= xmin_SidePlane)                 {Double_t a_p[3] = {xmin_SidePlane, y_SidePlane, z};}
      else if(x_th1 <= x_TopPlane)                     {Double_t a_p[3] = {x_TopPlane, y_TopPlane, z};}
      else if(z >= z_PCB)                              {Double_t a_p[3] = {x_th1, y_th1, z_PCB};}
      else if(z <= -z_PCB)                             {Double_t a_p[3] = {x_th1, y_th1, -z_PCB};}
      else                                             {Double_t a_p[3] = {x_th1, y_th1, z};}
   }else if(min_num == 6){
      if     (x_th2 >= xmin_SidePlane && z >= z_PCB)   {Double_t a_p[3] = {xmin_SidePlane, -y_SidePlane, z_PCB};}
      else if(x_th2 >= xmin_SidePlane && z <= -z_PCB)  {Double_t a_p[3] = {xmin_SidePlane, -y_SidePlane, -z_PCB};}
      else if(x_th2 <= x_TopPlane && z >= z_PCB)       {Double_t a_p[3] = {x_TopPlane, -y_TopPlane, z_PCB};}
      else if(x_th2 <= x_TopPlane && z <= -z_PCB)      {Double_t a_p[3] = {x_TopPlane, -y_TopPlane, -z_PCB};}
      else if(x_th2 >= xmin_SidePlane)                 {Double_t a_p[3] = {xmin_SidePlane, -y_SidePlane, z};}
      else if(x_th2 <= x_TopPlane)                     {Double_t a_p[3] = {x_TopPlane, -y_TopPlane, z};}
      else if(z >= z_PCB)                              {Double_t a_p[3] = {x_th2, y_th2, z_PCB};}
      else if(z <= -z_PCB)                             {Double_t a_p[3] = {x_th2, y_th2, -z_PCB};}
      else                                             {Double_t a_p[3] = {x_th2, y_th2, z};}
   }else if(min_num == 7){
      if     (x >= xmin_SidePlane && y >= y_SidePlane) {Double_t a_p[3] = {xmin_SidePlane, y_SidePlane, z_PCB};}
      else if(x >= xmin_SidePlane && y <= -y_SidePlane){Double_t a_p[3] = {xmin_SidePlane, -y_SidePlane, z_PCB};}
      else if(x <= x_TopPlane && y >= y_TopPlane)      {Double_t a_p[3] = {x_TopPlane, y_TopPlane, z_PCB};}
      else if(x <= x_TopPlane && y <= -y_TopPlane)     {Double_t a_p[3] = {x_TopPlane, -y_TopPlane, z_PCB};}
      else if(x >= xmin_SidePlane)                     {Double_t a_p[3] = {xmin_SidePlane, y, z_PCB};}
      else if(x <= x_TopPlane)                         {Double_t a_p[3] = {x_TopPlane, y, z_PCB};}
      else if(Abs(y) <= y_th)                          {Double_t a_p[3] = {x, y, z_PCB};}
      else if(y_th1 >= y_SidePlane)                    {Double_t a_p[3] = {xmin_SidePlane, y_SidePlane, z_PCB};}
      else if(y_th2 <= -y_SidePlane)                   {Double_t a_p[3] = {xmin_SidePlane, -y_SidePlane, z_PCB};}
      else if(y >= y_th)                               {Double_t a_p[3] = {x_th1, y_th1, z_PCB};}
      else if(y <= -y_th)                              {Double_t a_p[3] = {x_th2, y_th2, z_PCB};}
   }else if(min_num == 8){
      if     (x >= xmin_SidePlane && y >= y_SidePlane) {Double_t a_p[3] = {xmin_SidePlane, y_SidePlane, -z_PCB};}
      else if(x >= xmin_SidePlane && y <= -y_SidePlane){Double_t a_p[3] = {xmin_SidePlane, -y_SidePlane, -z_PCB};}
      else if(x <= x_TopPlane && y >= y_TopPlane)      {Double_t a_p[3] = {x_TopPlane, y_TopPlane, -z_PCB};}
      else if(x <= x_TopPlane && y <= -y_TopPlane)     {Double_t a_p[3] = {x_TopPlane, -y_TopPlane, -z_PCB};}
      else if(x >= xmin_SidePlane)                     {Double_t a_p[3] = {xmin_SidePlane, y, -z_PCB};}
      else if(x <= x_TopPlane)                         {Double_t a_p[3] = {x_TopPlane, y, -z_PCB};}
      else if(Abs(y) <= y_th)                          {Double_t a_p[3] = {x, y, -z_PCB};}
      else if(y_th1 >= y_SidePlane)                    {Double_t a_p[3] = {xmin_SidePlane, y_SidePlane, -z_PCB};}
      else if(y_th2 <= -y_SidePlane)                   {Double_t a_p[3] = {xmin_SidePlane, -y_SidePlane, -z_PCB};}
      else if(y >= y_th)                               {Double_t a_p[3] = {x_th1, y_th1, -z_PCB};}
      else if(y <= -y_th)                              {Double_t a_p[3] = {x_th2, y_th2, -z_PCB};}
   }

   r = m;
   a_prime = TVector3(a_p[0],a_p[1],a_p[2]);
}

void CalcML(TVector3 a[], TVectorD q, TMatrixD R, TMatrixD W, TMatrixD &M, TMatrixD &L){
   Long64_t N = sizeof(a)/sizeof(TVector3);
   TVector3 a_prime[N];
   
   Double_t r;

   Double_t xi[3][4];
   Double_t v[3];

   for(Int_t alpha=0; alpha<N ; alpha++){
      a[alpha] *= R * a[alpha]; 
      Dis(a[alpha],a_prime[alpha],r);

      xi[0][0] =  a_prime[alpha].X() - a[alpha].X();
      xi[1][0] =  0.;
      xi[2][0] = -a_prime[alpha].Z() - a[alpha].Z();
      xi[3][0] =  a_prime[alpha].Y() + a[alpha].Y();

      xi[0][1] =  a_prime[alpha].Y() - a[alpha].Y();
      xi[1][1] =  a_prime[alpha].Z() + a[alpha].Z();
      xi[2][1] =  0.;
      xi[3][1] = -a_prime[alpha].X() - a[alpha].X();

      xi[0][2] =  a_prime[alpha].Z() - a[alpha].Z();
      xi[1][2] = -a_prime[alpha].Y() - a[alpha].Y();
      xi[2][2] =  a_prime[alpha].X() + a[alpha].X();
      xi[3][2] =  0.;

      for(Int_t i=0; i<4; i++){
         for(Int_t j=0; j<4; j++){
            for(Int_t k=0; k<3; k++){
               for(Int_t l=0; l<3; l++){
                  M(i,j) += xi[i][k] * xi[j][l] * W(k,l);
               }
            }
         }
      }
      for(Int_t k=0; k<3; k++){
         v[k] = 0.;
         for(Int_t l=0; l<3; l++){
            for(Int_t i=0; i<4; i++){
               v[k] += W(k,l) * xi[i][l] * q(i);
            }
         }
      }
      for(Int_t k=0; k<3; k++){
         for(Int_t l=0; l<3; l++){
            L += v[k] * v[l] * V0_xi[k][l];
         }
      }
   }
}
