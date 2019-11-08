#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TMath.h>

/*TTree *t;
Int_t n0;
Double_t x0,y0,z0,yCut;

void circle(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    Double_t chisq=0;
    Double_t xx;
    for(Int_t i=0; i<n0; i++){
       t->GetEntry(i);
       if(y0 < yCut){continue;}
       else{
          xx = x0 - par[0];
          chisq += xx * xx + y0 * y0 - par[1];
       }
    }
    f = chisq;
}*/

void kiridashi() {
   const Double_t xCenter = 1.84934/(2*0.00249557);
   const Double_t yCutMin = 7.;
   const Double_t yCutMax = 340.;
   const Double_t Slope = -0.192;
   const Double_t Theta = TMath::Pi()/2 - TMath::ATan(0.192);
   const Double_t Slope2 = 5.6;

   Double_t x,y,z;
   Double_t phi,z2,phi2;
   Double_t x0,y0,z0,phi0,z20,phi20;
   Double_t xd,yd,zd;

   TFile *fout = new TFile("kiridashi.root","recreate");
   Int_t j = 23;
   TTree *t2[23];
   TTree *t3[256];
   TTree *t0 = new TTree("t0","t0");
   t0->Branch("x0", &xd); t0->Branch("y0", &yd); t0->Branch("z0", &zd);
   TTree *t1 = new TTree("t1","t1");
   t1->Branch("x", &x0); t1->Branch("y", &y0);
   t1->Branch("z", &z0); t1->Branch("phi", &phi0);
   t1->Branch("z2", &z20); t1->Branch("phi2", &phi20);

   TFile *data = new TFile("SPXScanReadData_190318_US.root");
   TTree *td = (TTree*)data->Get("tree");
   td->SetBranchAddress("x",&xd); td->SetBranchAddress("y",&yd);
   td->SetBranchAddress("z",&zd);
   Int_t nd = td->GetEntries();

   for(Int_t i=0; i<nd; i++){
      td->GetEntry(i);
      t0->Fill();
      if(yd <= yCutMin){continue;}
      else if(yd > yCutMax){continue;}
      else{
         x0 = xd - xCenter;
         y0 = yd - yCutMin;
         z0 = zd;
         phi0 = -(TMath::ATan(x0/y0) * TMath::RadToDeg() - 90); // ATan(y/x)だとπ/2を跨いでしまう。
         z20 = z0 * TMath::Cos(Theta) + phi0 * TMath::Sin(Theta); // カウンターを縦にした。
         phi20 = -z0 * TMath::Sin(Theta) + phi0 * TMath::Cos(Theta);
         t1->Fill();
      }
   }

   fout->cd();
   //t0->Write();
   //t1->Write();
   for(Int_t i=0; i<23; i++){
      t2[i] = new TTree(Form("t2_%d",i),Form("t2_%d",i));
      t2[i]->Branch("z2",&z2);
      t2[i]->Branch("phi2",&phi2);
      t2[i]->Branch("x",&x);
      t2[i]->Branch("y",&y);
      t2[i]->Branch("z",&z0);
      t2[i]->Branch("phi", &phi0);
   }
   for(Int_t i=0; i<256; i++){
      t3[i] = new TTree(Form("t3_%d",i),Form("t3_%d",i));
      t3[i]->Branch("z2",&z2);
      t3[i]->Branch("phi2",&phi2);
      t3[i]->Branch("x",&x);
      t3[i]->Branch("y",&y);
      t3[i]->Branch("z",&z0);
      t3[i]->Branch("phi", &phi0);
   }

   Double_t Phi1,Phi2;
   Int_t n1 = t1->GetEntries();
   for(Int_t i=0; i<23; i++){
      for(Int_t I=0; I<n1; I++){
         t1->GetEntry(I);
         Phi1 = Slope2 * z20 - 900 + i * 114;
         Phi2 = Slope2 * z20 - 740 + i * 114;
         if(phi20 > Phi1 && phi20 < Phi2){
            z2 = z20;
            phi2 = phi20;
            x = x0 + xCenter;
            y = y0 + yCutMin;
            t2[i]->Fill();
         }else{continue;}
      }
   }
   for(Int_t i=0; i<23; i++){
      t2[i]->Write();
   }
   Double_t ZMin = 0.;
   Int_t Nt = 0 ;
   Nt = t2[0]->GetEntries();
   for(Int_t k=0; k<Nt; k++){
      t2[0]->GetEntry(k);
      if(z2 <= 162. && z2 >= 155.5 && x > 0.){t3[0]->Fill();}
      if(z2 <= 155.5 && z2 >= 152.)          {t3[1]->Fill();}
   }
   /* 2~5 */
   Nt = t2[1]->GetEntries();
   for(Int_t k=0; k<Nt; k++){
      t2[1]->GetEntry(k);
      if(z2 <= 150.5 && z2 >= 145.1)          {t3[16]->Fill();}
      if(z2 <= 145.1 && z2 >= 142. && x > 45.){t3[17]->Fill();}
      if(z2 <= 141. && z2 >= 137. && y > 105.){t3[2]->Fill();}
      if(z2 <= 136. && z2 >= 132. && y > 150.){t3[3]->Fill();}
   }
   /* 6~11 */
   Nt = t2[2]->GetEntries();
   for(Int_t k=0; k<Nt; k++){
      t2[2]->GetEntry(k);
      if(z2 >= 111. && z2 <= 115.2 && z0 > -80. && y < 320.){t3[5]->Fill();}
      if(z2 >= 116.5 && z2 <= 120.5 && z0 > -100.)          {t3[4]->Fill();}
      if(z2 >= 121.5 && z2 <= 125. && y > 160.)             {t3[19]->Fill();}
      if(z2 >= 127. && z2 <= 130. && y > 100.)              {t3[18]->Fill();}
      if(z2 >= 131.5 && z2 <= 135.2 && y < 150.)            {t3[33]->Fill();}
      if(z2 >= 135.5 && z2 <= 140. && y < 105.)             {t3[32]->Fill();}
   }
   /* 12~19 */
   Nt = t2[3]->GetEntries();
   for(Int_t k=0; k<Nt; k++){
      t2[3]->GetEntry(k);
      if(z2 >= 90. && z2 <= 94.5 && phi2 < 100.)             {t3[7]->Fill();}
      if(z2 >= 96. && z2 <= 100. && phi2 < 120.)             {t3[6]->Fill();}
      if(z2 >= 101. && z2 <= 105. && phi2 < 160. && x < 280.){t3[21]->Fill();}
      if(z2 >= 106. && z2 <= 110. && y > 196.)               {t3[20]->Fill();}
      if(z2 >= 111. && z2 <= 115. && z0 < -75. && y > 150.)  {t3[35]->Fill();}
      if(z2 >= 116. && z2 <= 120. && y > 100. && y < 195.)   {t3[34]->Fill();}
      if(z2 >= 121. && z2 <= 124.5 && y < 150.)              {t3[49]->Fill();}
      if(z2 >= 124.5 && z2 <= 130. && y < 100.)              {t3[48]->Fill();}
   }
   /* 20~29 */
   Nt = t2[4]->GetEntries();
   for(Int_t k=0; k<Nt; k++){
      t2[4]->GetEntry(k);
      if(z2 >= 68.6 && z2 <= 73.6 && z0 > -70.)                   {t3[9]->Fill();}
      if(z2 >= 74.6 && z2 <= 78.6 && phi2 < 120.)                 {t3[8]->Fill();}
      if(z2 >= 80. && z2 <= 83.8 && phi2 < 150.)                  {t3[23]->Fill();}
      if(z2 >= 85.5 && z2 <= 89. && phi2 < 175.)                  {t3[22]->Fill();}
      if(z2 >= 90. && z2 <= 94. && phi2 < 215. && phi2 > 100.)    {t3[37]->Fill();}
      if(z2 >= 95. && z2 <= 99. && phi2 < 235. && phi2 > 120.)    {t3[36]->Fill();}
      if(z2 >= 100.4 && z2 <= 103.4 && phi2 > 155. && phi2 < 270.){t3[51]->Fill();}
      if(z2 >= 105.5 && z2 <= 109. && phi2 < 300. && phi2 > 180.) {t3[50]->Fill();}
      if(z2 >= 111. && z2 <= 114. && phi2 > 210.)                 {t3[65]->Fill();}
      if(z2 >= 114.7 && z2 <= 119. && phi2 > 240.)                {t3[64]->Fill();}
   }
   /* 30~41 */
   Nt = t2[5]->GetEntries();
   for(Int_t k=0; k<Nt; k++){
      t2[5]->GetEntry(k);
      if(z2 >= 47.8 && z2 <= 52.2 && phi2 < 85.)                     {t3[11]->Fill();}
      if(z2 >= 53. && z2 <= 57.4 && phi2 < 115.)                     {t3[10]->Fill();}
      if(z2 >= 58. && z2 <= 63. && phi2 < 148.)                      {t3[25]->Fill();}
      if(z2 >= 64. && z2 <= 67.6 && phi2 < 175.)                     {t3[24]->Fill();}
      if(z2 >= 68.6 && z2 <= 73. && phi2 < 200. && phi2 > 90.)       {t3[39]->Fill();}
      if(z2 >= 74.6  && z2 <= 78. && phi2 < 230. && phi2 > 120.)     {t3[38]->Fill();}
      if(z2 >= 79. && z2 <= 83.4 && phi2 < 265. && phi2 > 150.)      {t3[53]->Fill();}
      if(z2 >= 84. && z2 <= 88. && phi2 < 290. && phi2 > 180.)       {t3[52]->Fill();}
      if(z2 >= 90.5 && z2 <= ZMin+93. && phi2 < 320. && phi2 > 210.) {t3[67]->Fill();}
      if(z2 >= 95.4 && z2 <= 98. && phi2 < 350. && phi2 > 230.)      {t3[66]->Fill();}
      if(z2 >= 100. && z2 <= 103.2 && x > 50. && y > 55. && phi2 > 265.){t3[81]->Fill();}
      if(z2 >= 103.2 && z2 <= 109. && phi2 > 295.)                   {t3[80]->Fill();}
   }
   /* 42~55 */
   Nt = t2[6]->GetEntries();
   ZMin = 26.;
   for(Int_t k=0; k<Nt; k++){
      t2[6]->GetEntry(k);
      if(z2 >= ZMin && z2 <= ZMin+4.)    {t3[13]->Fill();}
      if(z2 >= ZMin+6. && z2 <= ZMin+10.){t3[12]->Fill();}
      if(z2 >= ZMin+11. && z2 <= ZMin+15.){t3[27]->Fill();}
      if(z2 >= ZMin+16.5 && z2 <= ZMin+20.5){t3[26]->Fill();}
      if(z2 >= ZMin+22. && z2 <= ZMin+26.){t3[41]->Fill();}
      if(z2 >= ZMin+28. && z2 <= ZMin+31.){t3[40]->Fill();}
      if(z2 >= ZMin+33. && z2 <= ZMin+36.){t3[55]->Fill();}
      if(z2 >= ZMin+38. && z2 <= ZMin+42.){t3[54]->Fill();}
      if(z2 >= ZMin+43. && z2 <= ZMin+46.){t3[69]->Fill();}
      if(z2 >= ZMin+48. && z2 <= ZMin+52.){t3[68]->Fill();}
      if(z2 >= ZMin+54. && z2 <= ZMin+57.){t3[83]->Fill();}
      if(z2 >= ZMin+59. && z2 <= ZMin+62.){t3[82]->Fill();}
      if(z2 >= 89.6 && z2 <= 92.6 && phi2 > 320. && phi2 < 430.){t3[97]->Fill();}
      if(z2 >= 93. && z2 <= 98. && phi2 > 350. && y < 105.)  {t3[96]->Fill();}
   }
   /* 56~71 */
   Nt = t2[7]->GetEntries();
   ZMin = 6.;
   for(Int_t k=0; k<Nt; k++){
      t2[7]->GetEntry(k);
      if(z2 >= 5.7 && z2 <= 10.7 && phi2 < 70.) {t3[15]->Fill();}
      if(z2 >= 11. && z2 <= 15.5 && phi2 < 100.){t3[14]->Fill();}
      if(z2 >= ZMin+10. && z2 <= ZMin+14.){t3[29]->Fill();}
      if(z2 >= ZMin+16. && z2 <= ZMin+20.){t3[28]->Fill();}
      if(z2 >= ZMin+21. && z2 <= ZMin+24.){t3[43]->Fill();}
      if(z2 >= ZMin+26. && z2 <= ZMin+29.){t3[42]->Fill();}
      if(z2 >= ZMin+32. && z2 <= ZMin+35.){t3[57]->Fill();}
      if(z2 >= ZMin+37. && z2 <= ZMin+40.){t3[56]->Fill();}
      if(z2 >= ZMin+42. && z2 <= ZMin+46.){t3[71]->Fill();}
      if(z2 >= ZMin+48. && z2 <= ZMin+52.){t3[70]->Fill();}
      if(z2 >= ZMin+53. && z2 <= ZMin+56.){t3[85]->Fill();}
      if(z2 >= ZMin+58. && z2 <= ZMin+61.){t3[84]->Fill();}
      if(z2 >= ZMin+62. && z2 <= ZMin+66.){t3[99]->Fill();}
      if(z2 >= ZMin+68. && z2 <= ZMin+72.){t3[98]->Fill();}
      if(z2 >= 78.5 && z2 <= 82.3 && phi2 < 485. && phi2 > 380.){t3[113]->Fill();}
      if(z2 >= 82.4 && z2 <= 88. && phi2 > 405. && x >0.)                {t3[112]->Fill();}
   }
   /* 72~87 */
   Nt = t2[8]->GetEntries();
   ZMin = -4.;
   for(Int_t k=0; k<Nt; k++){
      t2[8]->GetEntry(k);
      if(z2 >= -4.2 && z2 <= -0.2 && phi2 < 122.)     {t3[31]->Fill();}
      if(z2 >= 0.5 && z2 <= 4. && phi2 < 160.)        {t3[30]->Fill();}
      if(z2 >= ZMin+10. && z2 <= ZMin+14.){t3[45]->Fill();}
      if(z2 >= ZMin+15. && z2 <= ZMin+19.){t3[44]->Fill();}
      if(z2 >= ZMin+20. && z2 <= ZMin+24.){t3[59]->Fill();}
      if(z2 >= ZMin+25. && z2 <= ZMin+29.){t3[58]->Fill();}
      if(z2 >= ZMin+31. && z2 <= ZMin+35.){t3[73]->Fill();}
      if(z2 >= ZMin+36. && z2 <= ZMin+40.){t3[72]->Fill();}
      if(z2 >= ZMin+41. && z2 <= ZMin+45.){t3[87]->Fill();}
      if(z2 >= ZMin+47. && z2 <= ZMin+50.){t3[86]->Fill();}
      if(z2 >= ZMin+52. && z2 <= ZMin+56.){t3[101]->Fill();}
      if(z2 >= ZMin+57. && z2 <= ZMin+61.){t3[100]->Fill();}
      if(z2 >= ZMin+62. && z2 <= ZMin+66.){t3[115]->Fill();}
      if(z2 >= ZMin+68. && z2 <= ZMin+72.){t3[114]->Fill();}
      if(z2 >= 68.5 && z2 <= 72. && phi2 > 430. && phi2 < 560. && y > 55.){t3[129]->Fill();}
      if(z2 >= 72.1 && z2 <= 77.3 && phi2 > 460.){t3[128]->Fill();}
   }
   /* 88~103 */
   Nt = t2[9]->GetEntries();
   ZMin = -15.;
   for(Int_t k=0; k<Nt; k++){
      t2[9]->GetEntry(k);
      if(z2 >= -15. && z2 <= -10.5 && phi2 < 180.){t3[47]->Fill();}
      if(z2 >= -10. && z2 <= -6.5 && phi2 < 210.){t3[46]->Fill();}
      if(z2 >= ZMin+10. && z2 <= ZMin+14.){t3[61]->Fill();}
      if(z2 >= ZMin+15. && z2 <= ZMin+19.){t3[60]->Fill();}
      if(z2 >= ZMin+20. && z2 <= ZMin+24.){t3[75]->Fill();}
      if(z2 >= ZMin+26. && z2 <= ZMin+30.){t3[74]->Fill();}
      if(z2 >= ZMin+31. && z2 <= ZMin+35.){t3[89]->Fill();}
      if(z2 >= ZMin+36. && z2 <= ZMin+40.){t3[88]->Fill();}
      if(z2 >= ZMin+42. && z2 <= ZMin+46.){t3[103]->Fill();}
      if(z2 >= ZMin+47. && z2 <= ZMin+51.){t3[102]->Fill();}
      if(z2 >= ZMin+52. && z2 <= ZMin+56.){t3[117]->Fill();}
      if(z2 >= ZMin+57. && z2 <= ZMin+61.){t3[116]->Fill();}
      if(z2 >= ZMin+63. && z2 <= ZMin+67.){t3[131]->Fill();}
      if(z2 >= ZMin+69. && z2 <= ZMin+72.){t3[130]->Fill();}
      if(z2 >= 58. && z2 <= 61.6 && phi2 > 480. && phi2 <595. && x >60.){t3[145]->Fill();}
      if(z2 >= 62. && z2 <= 67. && phi2 > 510.){t3[144]->Fill();}
   }
   /* 104~119 */
   Nt = t2[10]->GetEntries();
   ZMin = -25.;
   for(Int_t k=0; k<Nt; k++){
      t2[10]->GetEntry(k);
      if(z2 >= -25. && z2 <= -20.5 && phi2 < 230.){t3[63]->Fill();}
      if(z2 >= -20.3 && z2 <= -17. && phi2 < 260.){t3[62]->Fill();}
      if(z2 >= ZMin+10. && z2 <= ZMin+14.){t3[77]->Fill();}
      if(z2 >= ZMin+15. && z2 <= ZMin+19.){t3[76]->Fill();}
      if(z2 >= ZMin+20. && z2 <= ZMin+24.){t3[91]->Fill();}
      if(z2 >= ZMin+26. && z2 <= ZMin+30.){t3[90]->Fill();}
      if(z2 >= ZMin+31. && z2 <= ZMin+35.){t3[105]->Fill();}
      if(z2 >= ZMin+36. && z2 <= ZMin+40.){t3[104]->Fill();}
      if(z2 >= ZMin+41. && z2 <= ZMin+45.){t3[119]->Fill();}
      if(z2 >= ZMin+47. && z2 <= ZMin+51.){t3[118]->Fill();}
      if(z2 >= ZMin+52. && z2 <= ZMin+56.){t3[133]->Fill();}
      if(z2 >= ZMin+57. && z2 <= ZMin+61.){t3[132]->Fill();}
      if(z2 >= ZMin+63. && z2 <= ZMin+67.){t3[147]->Fill();}
      if(z2 >= ZMin+68. && z2 <= ZMin+72.){t3[146]->Fill();}
      if(z2 >= 48. && z2 <= 51.2 && phi2 > 540.){t3[161]->Fill();}
      if(z2 >= 52. && z2 <= 56.5 && phi2 > 570.){t3[160]->Fill();}
   }
   /* 120~135 */
   Nt = t2[11]->GetEntries();
   ZMin = -36.;
   for(Int_t k=0; k<Nt; k++){
      t2[11]->GetEntry(k);
      if(z2 >= -36. && z2 <= -31. && phi2 < 290.)  {t3[79]->Fill();}
      if(z2 >= -30.7 && z2 <= -26.7 && phi2 < 315.){t3[78]->Fill();}
      if(z2 >= ZMin+10. && z2 <= ZMin+14.){t3[93]->Fill();}
      if(z2 >= ZMin+16. && z2 <= ZMin+20.){t3[92]->Fill();}
      if(z2 >= ZMin+21. && z2 <= ZMin+25.){t3[107]->Fill();}
      if(z2 >= ZMin+26. && z2 <= ZMin+29.){t3[106]->Fill();}
      if(z2 >= ZMin+31. && z2 <= ZMin+35.){t3[121]->Fill();}
      if(z2 >= ZMin+37. && z2 <= ZMin+41.){t3[120]->Fill();}
      if(z2 >= ZMin+42. && z2 <= ZMin+46.){t3[135]->Fill();}
      if(z2 >= ZMin+48. && z2 <= ZMin+52.){t3[134]->Fill();}
      if(z2 >= ZMin+53. && z2 <= ZMin+57.){t3[149]->Fill();}
      if(z2 >= ZMin+58. && z2 <= ZMin+62.){t3[148]->Fill();}
      if(z2 >= ZMin+64. && z2 <= ZMin+67.){t3[163]->Fill();}
      if(z2 >= ZMin+68. && z2 <= ZMin+72.){t3[162]->Fill();}
      if(z2 >= 37.6 && z2 <= 40.6 && phi2 > 595.&& y > 55.){t3[177]->Fill();}
      if(z2 >= 41.2 && z2 <= 46. && phi2 > 620. && phi2 < 720.) {t3[176]->Fill();}
   }
   /* 136~151 */
   Nt = t2[12]->GetEntries();
   ZMin = -46.;
   for(Int_t k=0; k<Nt; k++){
      t2[12]->GetEntry(k);
      if(z2 >= -46. && z2 <= -42. && phi2 < 340.){t3[95]->Fill();}
      if(z2 >= -41.3 && z2 <= -37.5 && phi2 < 375.){t3[94]->Fill();}
      if(z2 >= ZMin+10. && z2 <= ZMin+14.){t3[109]->Fill();}
      if(z2 >= ZMin+15. && z2 <= ZMin+19.){t3[108]->Fill();}
      if(z2 >= ZMin+20. && z2 <= ZMin+24.){t3[123]->Fill();}
      if(z2 >= ZMin+26. && z2 <= ZMin+30.){t3[122]->Fill();}
      if(z2 >= ZMin+31. && z2 <= ZMin+35.){t3[137]->Fill();}
      if(z2 >= ZMin+36. && z2 <= ZMin+40.){t3[136]->Fill();}
      if(z2 >= ZMin+41. && z2 <= ZMin+45.){t3[151]->Fill();}
      if(z2 >= ZMin+47. && z2 <= ZMin+51.){t3[150]->Fill();}
      if(z2 >= ZMin+52. && z2 <= ZMin+56.){t3[165]->Fill();}
      if(z2 >= ZMin+58. && z2 <= ZMin+61.){t3[164]->Fill();}
      if(z2 >= ZMin+62. && z2 <= ZMin+66.){t3[179]->Fill();}
      if(z2 >= ZMin+68. && z2 <= ZMin+72.){t3[178]->Fill();}
      if(z2 >= 27.9 && z2 <= 30.2 && phi2 > 645.){t3[193]->Fill();}
      if(z2 >= 31.2 && z2 <= 35.5 && phi2 > 675.){t3[192]->Fill();}
   }
   /* 152~167 */
   Nt = t2[13]->GetEntries();
   ZMin = -58.;
   for(Int_t k=0; k<Nt; k++){
      t2[13]->GetEntry(k);
      if(z2 >= -56.8 && z2 <= -51.7 && phi2 < 400.){t3[111]->Fill();}
      if(z2 >= -51.6 && z2 <= -48.5 && phi2 < 425.){t3[110]->Fill();}
      if(z2 >= ZMin+11. && z2 <= ZMin+15.){t3[125]->Fill();}
      if(z2 >= ZMin+16. && z2 <= ZMin+20.){t3[124]->Fill();}
      if(z2 >= ZMin+22. && z2 <= ZMin+25.){t3[139]->Fill();}
      if(z2 >= ZMin+26. && z2 <= ZMin+30.){t3[138]->Fill();}
      if(z2 >= ZMin+32. && z2 <= ZMin+36.){t3[153]->Fill();}
      if(z2 >= ZMin+38. && z2 <= ZMin+42.){t3[152]->Fill();}
      if(z2 >= ZMin+43. && z2 <= ZMin+47.){t3[167]->Fill();}
      if(z2 >= ZMin+48. && z2 <= ZMin+52.){t3[166]->Fill();}
      if(z2 >= ZMin+54. && z2 <= ZMin+58.){t3[181]->Fill();}
      if(z2 >= ZMin+59. && z2 <= ZMin+63.){t3[180]->Fill();}
      if(z2 >= ZMin+64. && z2 <= ZMin+68.){t3[195]->Fill();}
      if(z2 >= ZMin+69. && z2 <= ZMin+73.){t3[194]->Fill();}
      if(z2 >= 17. && z2 <= 20.2 && phi2 < 810. && phi2 > 700.){t3[209]->Fill();}
      if(z2 >= 23.2 && z2 <= 25. && phi2 < 834. && phi2 > 730. && x > 30.){t3[208]->Fill();}
   }
   /* 168~183 */
   Nt = t2[14]->GetEntries();
   ZMin = -68.;
   for(Int_t k=0; k<Nt; k++){
      t2[14]->GetEntry(k);
      if(z2 >= -67.2 && z2 <= -62.5 && phi2 < 450.)    {t3[127]->Fill();}
      if(z2 >= -62. && z2 <= -59.2 && phi2 < 475.)     {t3[126]->Fill();}
      if(z2 >= ZMin+11. && z2 <= ZMin+15.){t3[141]->Fill();}
      if(z2 >= ZMin+16. && z2 <= ZMin+20.){t3[140]->Fill();}
      if(z2 >= ZMin+21. && z2 <= ZMin+25.){t3[155]->Fill();}
      if(z2 >= ZMin+26. && z2 <= ZMin+30.){t3[154]->Fill();}
      if(z2 >= ZMin+32. && z2 <= ZMin+35.){t3[169]->Fill();}
      if(z2 >= ZMin+36. && z2 <= ZMin+40.){t3[168]->Fill();}
      if(z2 >= ZMin+42. && z2 <= ZMin+46.){t3[183]->Fill();}
      if(z2 >= ZMin+48. && z2 <= ZMin+52.){t3[182]->Fill();}
      if(z2 >= ZMin+53. && z2 <= ZMin+57.){t3[197]->Fill();}
      if(z2 >= ZMin+58. && z2 <= ZMin+62.){t3[196]->Fill();}
      if(z2 >= ZMin+63. && z2 <= ZMin+67.){t3[211]->Fill();}
      if(z2 >= ZMin+69. && z2 <= ZMin+73.){t3[210]->Fill();}
      if(z2 >= 6. && z2 <= 10. && phi2 > 755. && phi2 < 860.) {t3[225]->Fill();}
      if(z2 >= 11. && z2 <= 15. && phi2 > 785. && phi2 < 890. && x > 40.){t3[224]->Fill();}
   }
   /* 184~199 */
   Nt = t2[15]->GetEntries();
   ZMin = -79.;
   for(Int_t k=0; k<Nt; k++){
      t2[15]->GetEntry(k);
      if(z2 >= -77.2 && z2 <= -73. && phi2 < 500.){t3[143]->Fill();}
      if(z2 >= -72.4 && z2 <= -69.8 && phi2 < 540.){t3[142]->Fill();}
      if(z2 >= ZMin+11. && z2 <= ZMin+15.){t3[157]->Fill();}
      if(z2 >= ZMin+16. && z2 <= ZMin+20.){t3[156]->Fill();}
      if(z2 >= ZMin+21. && z2 <= ZMin+25.){t3[171]->Fill();}
      if(z2 >= ZMin+26. && z2 <= ZMin+30.){t3[170]->Fill();}
      if(z2 >= ZMin+32. && z2 <= ZMin+36.){t3[185]->Fill();}
      if(z2 >= ZMin+37. && z2 <= ZMin+41.){t3[184]->Fill();}
      if(z2 >= ZMin+43. && z2 <= ZMin+47.){t3[199]->Fill();}
      if(z2 >= ZMin+48. && z2 <= ZMin+52.){t3[198]->Fill();}
      if(z2 >= ZMin+53. && z2 <= ZMin+57.){t3[213]->Fill();}
      if(z2 >= ZMin+59. && z2 <= ZMin+63.){t3[212]->Fill();}
      if(z2 >= ZMin+64. && z2 <= ZMin+68.){t3[227]->Fill();}
      if(z2 >= ZMin+69. && z2 <= ZMin+73.){t3[226]->Fill();}
      if(z2 >= -4.5 && z2 <= -0.5 && phi2 > 810. && phi2 < 915. && x > 50.){t3[241]->Fill();}
      if(z2 >= 0.5 && z2 <= 4.3 && phi2 < 943. && phi2 > 840. && x > 40.)  {t3[240]->Fill();}
   }
   /* 200~213 */
   Nt = t2[16]->GetEntries();
   ZMin = -88.;
   for(Int_t k=0; k<Nt; k++){
      t2[16]->GetEntry(k);
      if(z2 >= -88. && z2 <= -83.5 && phi2 < 560.){t3[159]->Fill();}
      if(z2 >= -83. && z2 <= -80.5 && phi2 < 590.){t3[158]->Fill();}
      if(z2 >= ZMin+10. && z2 <= ZMin+14.){t3[173]->Fill();}
      if(z2 >= ZMin+15. && z2 <= ZMin+19.){t3[172]->Fill();}
      if(z2 >= ZMin+20. && z2 <= ZMin+24.){t3[187]->Fill();}
      if(z2 >= ZMin+25. && z2 <= ZMin+29.){t3[186]->Fill();}
      if(z2 >= ZMin+31. && z2 <= ZMin+35.){t3[201]->Fill();}
      if(z2 >= ZMin+36. && z2 <= ZMin+40.){t3[200]->Fill();}
      if(z2 >= ZMin+41. && z2 <= ZMin+45.){t3[215]->Fill();}
      if(z2 >= ZMin+47. && z2 <= ZMin+51.){t3[214]->Fill();}
      if(z2 >= ZMin+52. && z2 <= ZMin+56.){t3[229]->Fill();}
      if(z2 >= ZMin+58. && z2 <= ZMin+62.){t3[228]->Fill();}
      if(z2 >= ZMin+62.5 && z2 <= ZMin+66.5){t3[243]->Fill();}
      if(z2 >= ZMin+68. && z2 <= ZMin+72.){t3[242]->Fill();}
   }
   /* 214~225 */
   Nt = t2[17]->GetEntries();
   ZMin = -99.;
   for(Int_t k=0; k<Nt; k++){
      t2[17]->GetEntry(k);
      if(z2 >= -98.5 && z2 <= -94.5 && phi2 < 610.){t3[175]->Fill();}
      if(z2 >= -93.5 && z2 <= -90.5 && phi2 < 640.){t3[174]->Fill();}
      if(z2 >= ZMin+11. && z2 <= ZMin+14.5){t3[189]->Fill();}
      if(z2 >= ZMin+15.5 && z2 <= ZMin+19.){t3[188]->Fill();}
      if(z2 >= ZMin+20. && z2 <= ZMin+24.){t3[203]->Fill();}
      if(z2 >= ZMin+26. && z2 <= ZMin+30.){t3[202]->Fill();}
      if(z2 >= ZMin+31. && z2 <= ZMin+35.){t3[217]->Fill();}
      if(z2 >= ZMin+37. && z2 <= ZMin+41.){t3[216]->Fill();}
      if(z2 >= ZMin+43. && z2 <= ZMin+47.){t3[231]->Fill();}
      if(z2 >= ZMin+47.5 && z2 <= ZMin+51.){t3[230]->Fill();}
      if(z2 >= ZMin+52. && z2 <= ZMin+56.){t3[245]->Fill();}
      if(z2 >= ZMin+58. && z2 <= ZMin+62.){t3[244]->Fill();}
   }
   /* 226~235 */
   Nt = t2[18]->GetEntries();
   ZMin = -109.;
   for(Int_t k=0; k<Nt; k++){
      t2[18]->GetEntry(k);
      if(z2 >= -109. && z2 <= -105. && phi2 < 670.){t3[191]->Fill();}
      if(z2 >= -104. && z2 <= -101. && phi2 < 700.){t3[190]->Fill();}
      if(z2 >= ZMin+10. && z2 <= ZMin+14.){t3[205]->Fill();}
      if(z2 >= ZMin+15. && z2 <= ZMin+19.){t3[204]->Fill();}
      if(z2 >= ZMin+20. && z2 <= ZMin+24.){t3[219]->Fill();}
      if(z2 >= ZMin+25. && z2 <= ZMin+29.){t3[218]->Fill();}
      if(z2 >= ZMin+31. && z2 <= ZMin+35.){t3[233]->Fill();}
      if(z2 >= ZMin+36. && z2 <= ZMin+40.){t3[232]->Fill();}
      if(z2 >= ZMin+41. && z2 <= ZMin+45.){t3[247]->Fill();}
      if(z2 >= ZMin+47.5 && z2 <= ZMin+51.){t3[246]->Fill();}
   }
   /* 236~243 */
   Nt = t2[19]->GetEntries();
   ZMin = -119.5;
   for(Int_t k=0; k<Nt; k++){
      t2[19]->GetEntry(k);
      if(z2 >= -119.2 && z2 <= -115.5 && phi2 < 720.){t3[207]->Fill();}
      if(z2 >= -114.5 && z2 <= -111.5 && phi2 < 750.){t3[206]->Fill();}
      if(z2 >= ZMin+10.5 && z2 <= ZMin+13.5){t3[221]->Fill();}
      if(z2 >= ZMin+15.5 && z2 <= ZMin+19.){t3[220]->Fill();}
      if(z2 >= ZMin+20.5 && z2 <= ZMin+24.){t3[235]->Fill();}
      if(z2 >= ZMin+25.5 && z2 <= ZMin+29.5){t3[234]->Fill();}
      if(z2 >= ZMin+31. && z2 <= ZMin+34.){t3[249]->Fill();}
      if(z2 >= ZMin+37.5 && z2 <= ZMin+41.5){t3[248]->Fill();}
   }
   /* 244~249 */
   Nt = t2[20]->GetEntries();
   ZMin = -129.5;
   for(Int_t k=0; k<Nt; k++){
      t2[20]->GetEntry(k);
      if(z2 >= -129.5 && z2 <= -126.5 && phi2 < 769.){t3[223]->Fill();}
      if(z2 >= -124.5 && z2 <= -121.5 && phi2 < 800.){t3[222]->Fill();}
      if(z2 >= ZMin+9.5 && z2 <= ZMin+13.){t3[237]->Fill();}
      if(z2 >= ZMin+15. && z2 <= ZMin+18.){t3[236]->Fill();}
      if(z2 >= ZMin+19.5 && z2 <= ZMin+22.5){t3[251]->Fill();}
      if(z2 >= ZMin+26. && z2 <= ZMin+28.5){t3[250]->Fill();}
   }
   /* 250~253 */
   Nt = t2[21]->GetEntries();
   ZMin = -140.;
   for(Int_t k=0; k<Nt; k++){
      t2[21]->GetEntry(k);
      if(z2 >= -139.8 && z2 <= -136.6 && phi2 < 830.){t3[239]->Fill();}
      if(z2 >= -135.3 && z2 <= -132. && phi2 < 852.) {t3[238]->Fill();}
      if(z2 >= ZMin+10. && z2 <= ZMin+12.5){t3[253]->Fill();}
      if(z2 >= ZMin+14.5 && z2 <= ZMin+16.5){t3[252]->Fill();}
   }
   /* 254~255 */
   Nt = t2[22]->GetEntries();
   for(Int_t k=0; k<Nt; k++){
      t2[22]->GetEntry(k);
      if(z2 >= -150.4 && z2 <= -147.9 && phi2 < 877.){t3[255]->Fill();}
      if(z2 >= -145.7 && z2 <= -143.6 && phi2 < 898.){t3[254]->Fill();}
   }

   fout->cd();

   for(Int_t i=0; i<256; i++){
      t3[i]->Write();
   }
}
