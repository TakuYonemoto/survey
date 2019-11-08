#include "TROOT.h":q
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include <sstream>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TVector3.h>
#include <TMath.h>

using namespace TMath;

/// fit-box size ///
const Double_t hx = 61.9;   //mm
const Double_t Hy0 = 56.35; //mm
const Double_t hz2 = 3.1;    //mm
const Double_t hz = 5.25;  //mm
const Double_t tanyz = (TMath::Tan(TMath::ASin((hz-hz2)/Hy0)));
const Double_t hy = (hz-hz2)/tanyz/2; //mm
const Double_t S1 = (2*hz2 + 2*hz) * hy / 2; //mm^2
const Double_t S2 = 2*hz2 * 2*hx;            //mm^2
const Double_t S3 = 2*hx * Hy0;            //mm^2
const Double_t S4 = 2*hx * 2*hz;             //mm^2
const Double_t S  = 2*S1 + 2*S3 + S2 + S4; 

const Double_t hX = 60.;  //mm
const Double_t hY = 20.;  //mm
const Double_t hZ = 3.;  //mm

/// error ///
Double_t sigma = 0.025;  //mm

/// config ///

Double_t x,y,z,phi,z2,phi2,a,b,c;
Double_t cx2, cy2, cz2, cx_er2, cy_er2, cz_er2;
Double_t ta, tb, tc, ta_er, tb_er, tc_er;
Double_t ta2, tb2, tc2, ta_er2, tb_er2, tc_er2;
Double_t Gx,Gy,Gz,Ga,Gb,Gc,C1,C2,C3,C4;

void order()
{
   TFile *fo = new TFile("order.root","recreate");
   TTree *to = new TTree("tree","tree");
   to->Branch("dx",&x,"dx/D");
   to->Branch("dy",&y,"dy/D");
   to->Branch("dz",&z,"dz/D");
   to->Branch("da",&a,"da/D");
   to->Branch("db",&b,"db/D");
   to->Branch("dc",&c,"dc/D");

   TFile *f = new TFile("fit.root");
   TTree *tree = (TTree*)f->Get("tree");
   tree->SetBranchAddress("x",&x);
   tree->SetBranchAddress("y",&y);
   tree->SetBranchAddress("z",&z);
   tree->SetBranchAddress("a",&a);
   tree->SetBranchAddress("b",&b);
   tree->SetBranchAddress("c",&c);
   tree->GetEntry(1); to->Fill();//0
   tree->GetEntry(0); to->Fill();//1
   tree->GetEntry(3); to->Fill();//2
   tree->GetEntry(2); to->Fill();//3
   tree->GetEntry(7); to->Fill();//4
   tree->GetEntry(6); to->Fill();//5
   tree->GetEntry(13); to->Fill();//6
   tree->GetEntry(12); to->Fill();//7
   tree->GetEntry(21); to->Fill();//8
   tree->GetEntry(20); to->Fill();//9
   tree->GetEntry(31); to->Fill();//10
   tree->GetEntry(30); to->Fill();//11
   tree->GetEntry(43); to->Fill();//12
   tree->GetEntry(42); to->Fill();//13
   tree->GetEntry(57); to->Fill();//14
   tree->GetEntry(56); to->Fill();//15

   tree->GetEntry(5); to->Fill();//16
   tree->GetEntry(4); to->Fill();//17
   tree->GetEntry(9); to->Fill();//18
   tree->GetEntry(8); to->Fill();//19
   tree->GetEntry(15); to->Fill();//20
   tree->GetEntry(14); to->Fill();//21
   tree->GetEntry(23); to->Fill();//22
   tree->GetEntry(22); to->Fill();//23
   tree->GetEntry(33); to->Fill();//24
   tree->GetEntry(32); to->Fill();//25
   tree->GetEntry(45); to->Fill();//26
   tree->GetEntry(44); to->Fill();//27
   tree->GetEntry(59); to->Fill();//28
   tree->GetEntry(58); to->Fill();//29
   tree->GetEntry(73); to->Fill();//30
   tree->GetEntry(72); to->Fill();//31

   tree->GetEntry(11); to->Fill();//32
   tree->GetEntry(10); to->Fill();//33
   tree->GetEntry(17); to->Fill();//34
   tree->GetEntry(16); to->Fill();//35
   tree->GetEntry(25); to->Fill();//36
   tree->GetEntry(24); to->Fill();//37
   tree->GetEntry(35); to->Fill();//38
   tree->GetEntry(34); to->Fill();//39
   tree->GetEntry(47); to->Fill();//40
   tree->GetEntry(46); to->Fill();//41
   tree->GetEntry(61); to->Fill();//42
   tree->GetEntry(60); to->Fill();//43
   tree->GetEntry(75); to->Fill();//44
   tree->GetEntry(74); to->Fill();//45
   tree->GetEntry(89); to->Fill();//46
   tree->GetEntry(88); to->Fill();//47

   tree->GetEntry(19); to->Fill();//48
   tree->GetEntry(18); to->Fill();//49
   tree->GetEntry(27); to->Fill();//50
   tree->GetEntry(26); to->Fill();//51
   tree->GetEntry(37); to->Fill();//52
   tree->GetEntry(36); to->Fill();//53
   tree->GetEntry(49); to->Fill();//54
   tree->GetEntry(48); to->Fill();//55
   tree->GetEntry(63); to->Fill();//56
   tree->GetEntry(62); to->Fill();//57
   tree->GetEntry(77); to->Fill();//58
   tree->GetEntry(76); to->Fill();//59
   tree->GetEntry(91); to->Fill();//60
   tree->GetEntry(90); to->Fill();//61
   tree->GetEntry(105); to->Fill();//62
   tree->GetEntry(104); to->Fill();//63

   tree->GetEntry(29); to->Fill();//64
   tree->GetEntry(28); to->Fill();//65
   tree->GetEntry(39); to->Fill();//66
   tree->GetEntry(38); to->Fill();//67
   tree->GetEntry(51); to->Fill();//68
   tree->GetEntry(50); to->Fill();//69
   tree->GetEntry(65); to->Fill();//70
   tree->GetEntry(64); to->Fill();//71
   tree->GetEntry(79); to->Fill();//72
   tree->GetEntry(78); to->Fill();//73
   tree->GetEntry(93); to->Fill();//74
   tree->GetEntry(92); to->Fill();//75
   tree->GetEntry(107); to->Fill();//76
   tree->GetEntry(106); to->Fill();//77
   tree->GetEntry(121); to->Fill();//78
   tree->GetEntry(120); to->Fill();//79

   tree->GetEntry(41); to->Fill();//80
   tree->GetEntry(40); to->Fill();//81
   tree->GetEntry(53); to->Fill();//82
   tree->GetEntry(52); to->Fill();//83
   tree->GetEntry(67); to->Fill();//84
   tree->GetEntry(66); to->Fill();//85
   tree->GetEntry(81); to->Fill();//86
   tree->GetEntry(80); to->Fill();//87
   tree->GetEntry(95); to->Fill();//88
   tree->GetEntry(94); to->Fill();//89
   tree->GetEntry(109); to->Fill();//90
   tree->GetEntry(108); to->Fill();//91
   tree->GetEntry(123); to->Fill();//92
   tree->GetEntry(122); to->Fill();//93
   tree->GetEntry(137); to->Fill();//94
   tree->GetEntry(136); to->Fill();//95

   tree->GetEntry(55); to->Fill();//96
   tree->GetEntry(54); to->Fill();//97
   tree->GetEntry(69); to->Fill();//98
   tree->GetEntry(68); to->Fill();//99
   tree->GetEntry(83); to->Fill();//100
   tree->GetEntry(82); to->Fill();//101
   tree->GetEntry(97); to->Fill();//102
   tree->GetEntry(96); to->Fill();//103
   tree->GetEntry(111); to->Fill();//104
   tree->GetEntry(110); to->Fill();//105
   tree->GetEntry(125); to->Fill();//106
   tree->GetEntry(124); to->Fill();//107
   tree->GetEntry(139); to->Fill();//108
   tree->GetEntry(138); to->Fill();//109
   tree->GetEntry(153); to->Fill();//110
   tree->GetEntry(152); to->Fill();//111

   tree->GetEntry(71); to->Fill();//112
   tree->GetEntry(70); to->Fill();//113
   tree->GetEntry(85); to->Fill();//114
   tree->GetEntry(84); to->Fill();//115
   tree->GetEntry(99); to->Fill();//116
   tree->GetEntry(98); to->Fill();//117
   tree->GetEntry(113); to->Fill();//118
   tree->GetEntry(112); to->Fill();//119
   tree->GetEntry(127); to->Fill();//120
   tree->GetEntry(126); to->Fill();//121
   tree->GetEntry(141); to->Fill();//122
   tree->GetEntry(140); to->Fill();//123
   tree->GetEntry(155); to->Fill();//124
   tree->GetEntry(154); to->Fill();//125
   tree->GetEntry(169); to->Fill();//126
   tree->GetEntry(168); to->Fill();//127

   tree->GetEntry(87); to->Fill();//128
   tree->GetEntry(86); to->Fill();//129
   tree->GetEntry(101); to->Fill();//130
   tree->GetEntry(100); to->Fill();//131
   tree->GetEntry(115); to->Fill();//132
   tree->GetEntry(114); to->Fill();//133
   tree->GetEntry(129); to->Fill();//134
   tree->GetEntry(128); to->Fill();//135
   tree->GetEntry(143); to->Fill();//136
   tree->GetEntry(142); to->Fill();//137
   tree->GetEntry(157); to->Fill();//138
   tree->GetEntry(156); to->Fill();//139
   tree->GetEntry(171); to->Fill();//140
   tree->GetEntry(170); to->Fill();//141
   tree->GetEntry(185); to->Fill();//142
   tree->GetEntry(184); to->Fill();//143

   tree->GetEntry(103); to->Fill();//144
   tree->GetEntry(102); to->Fill();//145
   tree->GetEntry(117); to->Fill();//146
   tree->GetEntry(116); to->Fill();//147
   tree->GetEntry(131); to->Fill();//148
   tree->GetEntry(130); to->Fill();//149
   tree->GetEntry(145); to->Fill();//150
   tree->GetEntry(144); to->Fill();//151
   tree->GetEntry(159); to->Fill();//152
   tree->GetEntry(158); to->Fill();//153
   tree->GetEntry(173); to->Fill();//154
   tree->GetEntry(172); to->Fill();//155
   tree->GetEntry(187); to->Fill();//156
   tree->GetEntry(186); to->Fill();//157
   tree->GetEntry(201); to->Fill();//158
   tree->GetEntry(200); to->Fill();//159

   tree->GetEntry(119); to->Fill();//160
   tree->GetEntry(118); to->Fill();//161
   tree->GetEntry(133); to->Fill();//162
   tree->GetEntry(132); to->Fill();//163
   tree->GetEntry(147); to->Fill();//164
   tree->GetEntry(146); to->Fill();//165
   tree->GetEntry(161); to->Fill();//166
   tree->GetEntry(160); to->Fill();//167
   tree->GetEntry(175); to->Fill();//168
   tree->GetEntry(174); to->Fill();//169
   tree->GetEntry(189); to->Fill();//170
   tree->GetEntry(188); to->Fill();//171
   tree->GetEntry(203); to->Fill();//172
   tree->GetEntry(202); to->Fill();//173
   tree->GetEntry(215); to->Fill();//174
   tree->GetEntry(214); to->Fill();//175

   tree->GetEntry(135); to->Fill();//176
   tree->GetEntry(134); to->Fill();//177
   tree->GetEntry(149); to->Fill();//178
   tree->GetEntry(148); to->Fill();//179
   tree->GetEntry(163); to->Fill();//180
   tree->GetEntry(162); to->Fill();//181
   tree->GetEntry(177); to->Fill();//182
   tree->GetEntry(176); to->Fill();//183
   tree->GetEntry(191); to->Fill();//184
   tree->GetEntry(190); to->Fill();//185
   tree->GetEntry(205); to->Fill();//186
   tree->GetEntry(204); to->Fill();//187
   tree->GetEntry(217); to->Fill();//188
   tree->GetEntry(216); to->Fill();//189
   tree->GetEntry(227); to->Fill();//190
   tree->GetEntry(226); to->Fill();//191

   tree->GetEntry(151); to->Fill();//192
   tree->GetEntry(150); to->Fill();//193
   tree->GetEntry(165); to->Fill();//194
   tree->GetEntry(164); to->Fill();//195
   tree->GetEntry(179); to->Fill();//196
   tree->GetEntry(178); to->Fill();//197
   tree->GetEntry(193); to->Fill();//198
   tree->GetEntry(192); to->Fill();//199
   tree->GetEntry(207); to->Fill();//200
   tree->GetEntry(206); to->Fill();//201
   tree->GetEntry(219); to->Fill();//202
   tree->GetEntry(218); to->Fill();//203
   tree->GetEntry(228); to->Fill();//204
   tree->GetEntry(227); to->Fill();//205
   tree->GetEntry(237); to->Fill();//206
   tree->GetEntry(236); to->Fill();//207

   tree->GetEntry(167); to->Fill();//208
   tree->GetEntry(166); to->Fill();//209
   tree->GetEntry(181); to->Fill();//210
   tree->GetEntry(180); to->Fill();//211
   tree->GetEntry(195); to->Fill();//212
   tree->GetEntry(194); to->Fill();//213
   tree->GetEntry(209); to->Fill();//214
   tree->GetEntry(208); to->Fill();//215
   tree->GetEntry(221); to->Fill();//216
   tree->GetEntry(220); to->Fill();//217
   tree->GetEntry(231); to->Fill();//218
   tree->GetEntry(230); to->Fill();//219
   tree->GetEntry(239); to->Fill();//220
   tree->GetEntry(238); to->Fill();//221
   tree->GetEntry(245); to->Fill();//222
   tree->GetEntry(244); to->Fill();//223

   tree->GetEntry(183); to->Fill();//224
   tree->GetEntry(182); to->Fill();//225
   tree->GetEntry(197); to->Fill();//226
   tree->GetEntry(196); to->Fill();//227
   tree->GetEntry(211); to->Fill();//228
   tree->GetEntry(210); to->Fill();//229
   tree->GetEntry(223); to->Fill();//230
   tree->GetEntry(222); to->Fill();//231
   tree->GetEntry(233); to->Fill();//232
   tree->GetEntry(232); to->Fill();//233
   tree->GetEntry(241); to->Fill();//234
   tree->GetEntry(240); to->Fill();//235
   tree->GetEntry(247); to->Fill();//236
   tree->GetEntry(246); to->Fill();//237
   tree->GetEntry(251); to->Fill();//238
   tree->GetEntry(250); to->Fill();//239

   tree->GetEntry(199); to->Fill();//240
   tree->GetEntry(198); to->Fill();//241
   tree->GetEntry(213); to->Fill();//242
   tree->GetEntry(212); to->Fill();//243
   tree->GetEntry(225); to->Fill();//244
   tree->GetEntry(224); to->Fill();//245
   tree->GetEntry(235); to->Fill();//246
   tree->GetEntry(234); to->Fill();//247
   tree->GetEntry(243); to->Fill();//248
   tree->GetEntry(242); to->Fill();//249
   tree->GetEntry(249); to->Fill();//250
   tree->GetEntry(248); to->Fill();//251
   tree->GetEntry(253); to->Fill();//252
   tree->GetEntry(252); to->Fill();//253
   tree->GetEntry(255); to->Fill();//254
   tree->GetEntry(254); to->Fill();//255

   fo->cd();
   to->Write();
}
