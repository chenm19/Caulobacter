function param()
global p

p.clpxp=1;%um


%% Synthesis and degradation rate constants [units --> 1/min]
p.k1_pos=0.1;% %forward reaction: ClpXP+CpdR ->Complex1
p.k1_neg=0.1;% %backward reaction: Complex1->ClpXP+CpdR
p.ks_cpdr=0.6;%20;% 
p.kd_cpdr=1.2;%30;%
%%phosporylation
p.k2_pos=1;%6;%;%%%%sharpness of CpdRP%%%%
p.k2_neg=1;%6;%;%based on ratio of p/unp
%%
p.k3_pos=0.5;%forward reaction: complex1+RcdA ->Complex1
p.k3_neg=0.5;%backward reaction: complex1+RcdA
p.ks_rcda=0.023*20;%0.023*20;%0.023;
p.kd_rcda=0.017*20;%0.015*20;%0.017;%from Tyson lab
p.k5_pos=0.2;%
p.k5_neg=0.2;

%%
p.J1=4;%10%CpdR-Complex1

p.J3=2;%%0.5;%RcdA
p.J4=4;%0.4;%RcdA-Complex1



