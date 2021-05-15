function param()
global p

p.clpxp=1;%um

% update 3/13/2021 RX
%% Synthesis and degradation rate constants [units --> 1/min]
% p.k1_pos=0.1;% %forward reaction: ClpXP+CpdR ->Complex1
% p.k1_neg=0.1;% %backward reaction: Complex1->ClpXP+CpdR
% p.ks_cpdr=0.6;%20;% 
% p.kd_cpdr=1.2;%30;%
p.k1_pos=0.1*6;% %forward reaction: ClpXP+CpdR ->Complex1
p.k1_neg=3;%6;% %backward reaction: Complex1->ClpXP+CpdR         %2/18/2021RX
p.ks_cpdr=0.7;%;%20;% 
p.kd_cpdr=1.5;%30;%                                             %2/18/2021RX
%%phosporylation
p.dephos_CpdR=1;%1;%6;%;%%%%sharpness of CpdRP%%%%
p.phos_CpdR=1;%6;%;%based on ratio of p/unp                         %2/18/2021RX
 
p.JdCpdR=6;%10%CpdR-Complex1
p.naCpdRCtrA=2;%1;                                               %2/18/2021RX
p.JaCpdRCtrA =15;       %15                                         %2/18/2021RX
%%
p.k2_pos=1.1;%1%forward reaction: complex1+RcdA ->Complex1          %2/18/2021RX
p.k2_neg=1;%backward reaction: complex1+RcdA
p.ks_rcda=0.015*10;%%0.023;                                     %2/18/2021RX
p.kd_rcda=0.02*10;%%0.017;%from Tyson lab
p.k3_pos=140;%0.04;%%C2->C3
p.k3_neg=2;
%%
p.JaRcdACtrA=15;%%0.5;%RcdA
p.naRcdACtrA=2;
p.JdRcdA=2;%0.4;%RcdA-Complex1
%%cdG                                                       %2/18/2021RX
p.kscdG=0.01;p.kdcdG=1;
p.DGC=1500; p.PDE=7;           
p.JicdGcdG=0.2;%0.05; 
% p.JdcdG=0.2;
%%PleD
p.ksPleD=0.1; p.kdPleD=0.15;
p.naPleDCtrA=2; p.JaPleDCtrA=2.5;   %1;15                           %2/18/2021RX
p.phosPleD=2/50;   p.dephosPleD=2/50;   p.c=2;                                 %3/13/2021 RX
% p.ksPleD2=0;
%%PdeA
p.ksPdeA=0.02/2; p.kdPdeA=1/2;
p.naPdeACtrA=2;  p.JaPdeACtrA=5;
p.JdPdeA=5;