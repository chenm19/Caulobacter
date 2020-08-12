function dydt = ODEs5regCPLX(t, y)
%%variable map

global Ini Elong DNA Count hcori hCcrM hCtrA mCcrM mDnaA ...
    mGcrA mSciP mCtrA CcrM DnaA GcrA SciP CtrA Sup DivKp I II III tot;

global CPLX1 CpdR CpdRP CPLX2 RcdA CPLX3; 

global p;

%y(1) = [Ini]
%y(2) = [Elong]
%y(3) = [DNA]
%y(4) = [Count]
%y(5) = [hcori]
%y(6) = [hccrm]
%y(7) = [hctra]
%y(8) = [mCcrM]
%y(9) = [mDnaA]
%y(10) = [mGcrA]
%y(11) = [mSciP]
%y(12) = [mCtrA]
%y(13) = [CcrM]
%y(14) = [DnaA]
%y(15) = [GcrA]
%y(16) = [SciP]
%y(17) = [CtrA]
%y(18) = [Sup]
%y(19) = [DivKp]
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%&&&&&&&&%

% parameters for which cell we are tracking on
global CELLTYPE;

% cell type track
H = CELLTYPE; 
% parameters for the equations

% DNA replication, methylation
kelong = 0.95/160*1.1;
nelong = 4;
kmcori = 1.4;	Jmcori = 0.95;	nmcori = 4;
kmccrM = 1.4;	JmccrM = 0.95;	nmccrM = 4;
kmctrA = 1.4;	JmctrA = 0.95;	nmctrA = 4;

% DnaA, GcrA, CtrA, CcrM, SciP
ksCcrM = 0.1;  kdCcrM = 0.1;
ksDna = 0.1;   kdDna = 0.1;
ksGcrA = 0.1;  kdGcrA = 0.1;
ksSciP = 0.08;  kdSciP = 0.08;
ksCtrA = 0.024;  kdCtrA = 0.002;
kdCtrADivKp = 0.15; ndCtrADivKp = 2; JdCtrADivKp = 1;

ksI = 0.05; kdI = 0.05;
ksII = 0.05; kdII = 0.05;
ksIII = 0.1; kdIII = 0.1;
% DivKp

%mCcrM
ksmCcrM = 0.32;     kdmCcrM = 0.08;
JaCcrMCtrA = 3.5;   naCcrMCtrA = 2;
JiCcrMSciP = 3;   niCcrMSciP = 2;


%mDnaA
ksmDnaA = 0.055;   kdmDnaA = 0.015;
JiDnaAGcrA = 3.0;	niDnaAGcrA = 2;
JaDnaACtrA = 3.5;   naDnaACtrA = 2;

%mGcrA
ksmGcrA = 2.8;    kdmGcrA = 0.3;
JiGcrACtrA = 3.5;   niGcrACtrA = 2;
JaGcrADnaA = 1.5;   naGcrADnaA = 2;

%mSciP
ksmSciP = 0.5; kdmSciP = 0.06;
JaSciPCtrA = 3.5; naSciPCtrA = 2;

%mCtrA
ksmCtrA = 0.9;	kdmCtrA = 0.1; 

JaCtrACtrA = 3.5;	naCtrACtrA = 2;
JaCtrAGcrA = 3.0;	naCtrAGcrA = 2;
JiCtrACtrA = 3.5;	niCtrACtrA = 2;
JiCtrASciP = 3.0;	niCtrASciP = 2;

% end of parameters

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%&&&&&&&&%

% This section is used for simulations of mutants 
%if t > 120
%   kop = (0.0083 + 0.073)*1.5;
%	 ksCtrAP1 = 0;
%	 ksCtrAP2 = 0;
%end
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%&&&&&&&&%

% Differential equations of the model
dydt = zeros(tot, 1);  

% DNA replication
kaIni = 0.05;	JaIni = 1;		naIni = 4;

thetaCtrA = 1;	nthetaCtrA = 4;
thetaDnaA = 3;	nthetaDnaA = 4;

dydt(Ini) = kaIni*(JaIni^naIni/(JaIni^naIni + (y(CtrA)/thetaCtrA)^nthetaCtrA) ...
     * (y(DnaA)/thetaDnaA)^nthetaDnaA/(JaIni^naIni + (y(DnaA)/thetaDnaA)^nthetaDnaA) );

dydt(Elong) = kelong*y(Elong)^nelong/(y(Elong)^nelong + 0.05^nelong)*y(Count);

dydt(DNA) = kelong*y(Elong)^nelong/(y(Elong)^nelong + 0.05^nelong)*y(Count);

dydt(Count) = 0;

% DNA methylation
dydt(hcori) =  - kmcori*y(CcrM)^nmcori/(Jmcori^nmcori + y(CcrM)^nmcori)*y(hcori);
dydt(hCcrM) =  - kmccrM*y(CcrM)^nmctrA/(JmctrA^nmctrA + y(CcrM)^nmctrA)*y(hCcrM);
dydt(hCtrA) =  - kmctrA*y(CcrM)^nmccrM/(JmccrM^nmccrM + y(CcrM)^nmccrM)*y(hCtrA);

% mCcrM mDnaA mGcrA mSciP mCtrA
dydt(III) = ksIII*(y(CtrA)^naCcrMCtrA/(JaCcrMCtrA^naCcrMCtrA + y(CtrA)^naCcrMCtrA) ...
              * JiCcrMSciP^niCcrMSciP/(JiCcrMSciP^niCcrMSciP + y(SciP)^niCcrMSciP) ) ...
              * y(hCcrM) - kdIII*y(III);
dydt(mCcrM) = ksmCcrM*y(III) - kdmCcrM*y(mCcrM);
      
dydt(mDnaA) = ( ksmDnaA*JiDnaAGcrA^niDnaAGcrA/(JiDnaAGcrA^niDnaAGcrA + y(GcrA)^niDnaAGcrA) ...
              * y(CtrA)^naDnaACtrA/(JaDnaACtrA^naDnaACtrA + y(CtrA)^naDnaACtrA) ) ...
              * (2 - y(hcori)) - kdmDnaA*y(mDnaA);
          
dydt(mGcrA) = ( ksmGcrA*y(DnaA)^naGcrADnaA/(JaGcrADnaA^naGcrADnaA + y(DnaA)^naGcrADnaA) ...
              * JiGcrACtrA^niGcrACtrA/(JiGcrACtrA^niGcrACtrA + y(CtrA)^niGcrACtrA) ) ...
              - kdmGcrA*y(mGcrA);

dydt(mSciP) = ksmSciP*y(CtrA)^naSciPCtrA/(JaSciPCtrA^naSciPCtrA + y(CtrA)^naSciPCtrA) ...
              - kdmSciP*y(mSciP);

dydt(mCtrA) = (ksmCtrA*y(CtrA)^naCtrACtrA/(JaCtrACtrA^naCtrACtrA + y(CtrA)^naCtrACtrA) ...
              + ksmCtrA*y(GcrA)^naCtrAGcrA/(JaCtrAGcrA^naCtrAGcrA + y(GcrA)^naCtrAGcrA) ...
              *JiCtrACtrA^niCtrACtrA/(JiCtrACtrA^niCtrACtrA + y(CtrA)^niCtrACtrA) ...
              *JiCtrASciP^niCtrASciP/(JiCtrASciP^niCtrASciP + y(SciP)^niCtrASciP) ) ...
              *y(hCtrA) - kdmCtrA*y(mCtrA) ;

% CcrM DnaA GcrA SciP CtrA CtrA
dydt(CcrM) = ksCcrM*y(mCcrM) - kdCcrM*y(CcrM);
dydt(I) = ksI*y(mDnaA) - kdI*y(I);
dydt(DnaA) = ksDna*y(I) - kdDna*y(DnaA);

dydt(II) = ksII*y(mGcrA) - kdII*y(II);
dydt(GcrA) = ksGcrA*y(II) - kdGcrA*y(GcrA);
dydt(SciP) = ksSciP*y(mSciP) - kdSciP*y(SciP);


dydt(CtrA) = ksCtrA*y(mCtrA) - (kdCtrA + kdCtrADivKp*y(CPLX3)^ndCtrADivKp/...
             (JdCtrADivKp^ndCtrADivKp + y(CPLX3)^ndCtrADivKp) )*y(CtrA);

%          floor(t)

% sup
dydt(Sup) = 0.01;

% DivKp
dydt(DivKp) = 0;

%%
T=150;        % period of Caulobacter
t_d=rem(t,T); % return remainder after division t/T

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


%% cdG

cdG=0.4245*sin(pi/75*t_d+0.3976)+0.3872;
cdG(cdG<0)=0;
% Goodness of fit:
%   SSE: 0.1651
%   R-square: 0.8034
%   Adjusted R-square: 0.7248
%   RMSE: 0.1817
%%%????????
%%

dydt(CPLX1)=p.k1_pos*p.clpxp*y(CpdR)-p.k1_neg*y(CPLX1)...
    -p.k3_pos*y(CPLX1)*y(RcdA)+p.k3_neg*y(CPLX2);

dydt(CpdR)=p.ks_cpdr-p.kd_cpdr*y(CpdR)*y(CPLX1)/(y(CPLX1)+p.J1)...
    +p.k1_neg*y(CPLX1)-p.k1_pos*p.clpxp*y(CpdR)/(y(CpdR)+p.Km1)...
    +p.k2_pos*y(CpdRP)*p.divkp/(p.divkp+p.J2)-p.k2_neg*y(CpdR);

dydt(CpdRP)=p.k2_neg*y(CpdR)-p.k2_pos*y(CpdRP)*p.divkp/(p.divkp+p.J2);


dydt(CPLX2)=p.k3_pos*y(CPLX1)*y(RcdA)-p.k3_neg*y(CPLX2)+p.k5_neg*y(CPLX3)-p.k5_pos*cdG^2*y(CPLX2);
dydt(RcdA)=p.ks_rcda*y(RcdA)^1/(y(RcdA)^1+p.J3^1)-p.kd_rcda*y(RcdA)*y(CPLX1)/(y(CPLX1)+p.J4);


dydt(CPLX3)=-p.k5_neg*y(CPLX3)+p.k5_pos*cdG^2*y(CPLX2);

% end of equations
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%&&&&&&&&%

end
