function param()
global p

p.clpxp=1;%um


%% Synthesis and degradation rate constants [units --> 1/min]
p.k1_pos=0.2;%0.2; %forward reaction: ClpXP+CpdR ->Complex1
p.k1_neg=0.2;%20; %backward reaction: Complex1->ClpXP+CpdR
p.ks_cpdr=20;%10; %0.02;
p.kd_cpdr=30;%0.5;%0.01;
%%phosporylation
p.k2_pos=6;%2;%0.5;%%%%sharpness of CpdRP%%%%
p.k2_neg=6;%15;%0.8;%based on ratio of p/unp
%%
p.k3_pos=2;%forward reaction: complex1+RcdA ->Complex1
p.k3_neg=1;%backward reaction: complex1+RcdA
p.ks_rcda=0.023*20;%2.5;%0.023;
p.kd_rcda=0.015*20;%2;%0.017;%from Tyson lab
p.k4_pos=10;
p.k4_neg=5;
p.ks_cdg=0.5;
p.kd_cdg=1.5;
p.k5_pos=0.2;%0.2;%20;
p.k5_neg=0.2;

%%
p.J1=2;%10%CpdR-Complex1
p.J2=0.5;%0.5;
p.J3=3;%%0.5;%RcdA
p.J4=3;%0.4;%RcdA-Complex1
p.Km1=1;%5;
p.J5=0.5;
p.J6=2;
p.km2=0.1;
p.km3=0.06;
p.J7=0.3;
%%binding
p.kcpdr_b_f=0.1;
p.kcpdr_f_b=15;
p.kcpdrp_b_f=0.5;
p.kcpdrp_f_b=1;
p.kcdg_b_f=0.2;
p.kcdg_f_b=2;
