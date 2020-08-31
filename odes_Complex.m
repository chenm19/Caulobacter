%odes.m


function dydt=odes(t,Y)
global p;

Complex1=Y(1);
CpdR=Y(2);
CpdRP=Y(3);
Complex2=Y(4);
RcdA=Y(5);
Complex3=Y(6);
%%
T=150;%period of Caulobacter
t_d=rem(t,T); %return remainder after division t/T



%%
%% cdG
cdG= 0.3233 *sin(pi/75*t_d+0.3469 )+ 0.3363;
%complex1
dComplex1=p.k1_pos*p.clpxp*CpdR-p.k1_neg*Complex1...
    -p.k3_pos*Complex1*RcdA+p.k3_neg*Complex2;

%CpdR
dCpdR=p.ks_cpdr-p.kd_cpdr*CpdR*Complex1/(Complex1+p.J1)...
    +p.k1_neg*Complex1-p.k1_pos*p.clpxp*CpdR...
    +p.k2_pos*CpdRP-p.k2_neg*CpdR;
%CpdR~P
dCpdRP=p.k2_neg*CpdR-p.k2_pos*CpdRP;

%Complex2
dComplex2=p.k3_pos*Complex1*RcdA-p.k3_neg*Complex2+p.k5_neg*Complex3-p.k5_pos*cdG^2*Complex2;
%RcdA
dRcdA=p.ks_rcda*RcdA^1/(RcdA^1+p.J3^1)-p.kd_rcda*RcdA*Complex1/(Complex1+p.J4);
%Complex3
dComplex3=-p.k5_neg*Complex3+p.k5_pos*cdG^2*Complex2;


dydt=[dComplex1; dCpdR; dCpdRP; dComplex2; dRcdA; dComplex3];

