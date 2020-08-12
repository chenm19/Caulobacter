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

%DivKP
p.divkp=0.004218*t_d+ 0.5779;
p.divkp(p.divkp<0)=0;
%0.92
% Linear model Poly3:
  p.divkp = 1.64e-06*t_d.^3-0.0003057 *t_d.^2 + 0.02095*t_d +0.3282 ;
  
% Coefficients (with 95% confidence bounds):
%linear model poly1 divkp
p.divkp=0.004218*t_d+0.5779;
     
%RcdA
q1 =  -8.725e-08 ;
       q2 =   2.948e-05 ;
       q3 =   -0.003332  ;
       q4 =      0.1384 ;
       q5 =     -0.8569;
       p.rcda=q1*t_d.^4 + q2*t_d.^3 + q3*t_d.^2 + q4*t_d+ q5;



%%
%% cdG

% cdG=2.335*sin(0.02876*t_d-0.4084) + 1.955*sin(0.03501*t_d+2.189);%Goodness of fit:
% %   SSE: 0.03746
% %   R-square: 0.8904
% %   Adjusted R-square: 0.6163
% %   RMSE: 0.1369
%%%% make the last four points as 0.1
cdG=0.4245*sin(pi/75*t_d+0.3976)+0.3872;
cdG(cdG<0)=0;
% Goodness of fit:
%   SSE: 0.1651
%   R-square: 0.8034
%   Adjusted R-square: 0.7248
%   RMSE: 0.1817
%%%????????
%%
%xx=[20,40,60,80,100,120,140]
%yyy=[0.091,0.276,0.156,0.104,0.09,0.075,0.036]
% General model:
%      f(x) = a*sin(pi/75*x+b)+c
% Coefficients (with 95% confidence bounds):
%        a =     0.07759  (-0.008303, 0.1635)
%        b =     -0.4908  (-1.679, 0.6978)
%        c =      0.1167  (0.05373, 0.1797)
% 
% Goodness of fit:
%   SSE: 0.01424
%   R-square: 0.6118
%   Adjusted R-square: 0.4177
%   RMSE: 0.05967
% cdG= 0.07759 *sin(pi/75*t_d-0.4908 )+ 0.1167;
cdG= 0.3233 *sin(pi/75*t_d+0.3469 )+ 0.3363;
%[0    20    40    60    80   100   120]
%[0.2292    1.0000    0.5000    0.2833    0.2250    0.1625         0]
%%
%%%%%%%%%%Complex1 (ClpXP:CpdR)
% dComplex1=p.k1_pos*p.clpxp*CpdR/(CpdR+p.Km1)-p.k1_neg*Complex1;
% dComplex1=p.k1_pos*p.clpxp*CpdR-p.k1_neg*Complex1;
dComplex1=p.k1_pos*p.clpxp*CpdR-p.k1_neg*Complex1...
    -p.k3_pos*Complex1*RcdA+p.k3_neg*Complex2;
% dCpdR=p.ks_cpdr-p.kd_cpdr*CpdR*Complex1/(Complex1+p.J1)...
%     +p.k1_neg*Complex1-p.k1_pos*p.clpxp*CpdR/(CpdR+p.Km1)...
%     +p.k2_pos*CpdRP*p.divkp/(p.divkp+p.J2)-p.k2_neg*CpdR;
dCpdR=p.ks_cpdr-p.kd_cpdr*CpdR*Complex1/(Complex1+p.J1)...
    +p.k1_neg*Complex1-p.k1_pos*p.clpxp*CpdR...
    +p.k2_pos*CpdRP*p.divkp/(p.divkp+p.J2)-p.k2_neg*CpdR;
%     +p.k1_neg*Complex1-p.k1_pos*p.clpxp*CpdR/(CpdR+p.Km1)...
dCpdRP=p.k2_neg*CpdR-p.k2_pos*CpdRP*p.divkp/(p.divkp+p.J2);

% dydt=[dComplex1; dCpdR; dCpdRP];

% dComplex2=p.k3_pos*Complex1*RcdA-p.k3_neg*Complex2;
dComplex2=p.k3_pos*Complex1*RcdA-p.k3_neg*Complex2+p.k5_neg*Complex3-p.k5_pos*cdG^2*Complex2;
dRcdA=p.ks_rcda*RcdA^1/(RcdA^1+p.J3^1)-p.kd_rcda*RcdA*Complex1/(Complex1+p.J4);
% dRcdA=p.ks_rcda*RcdA^4/(RcdA^4+p.J3^4)-p.kd_rcda*RcdA*Complex1/(Complex1+p.J4);
% dRcdA=p.ks_rcda-p.kd_rcda*RcdA*Complex1/(Complex1+p.J4);

% dRcdA=p.ks_rcda*RcdA/(RcdA+p.J3)-p.kd_rcda*RcdA*Complex1;

% dydt=[dComplex1; dCpdR; dCpdRP; dComplex2; dRcdA];
%%

dComplex3=-p.k5_neg*Complex3+p.k5_pos*cdG^2*Complex2;
dydt=[dComplex1; dCpdR; dCpdRP; dComplex2; dRcdA; dComplex3];

