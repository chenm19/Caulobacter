function dydt = ODEs5regCPLX(t, y)
%%variable map

global Ini Elong DNA Count hcori hCcrM hCtrA mCcrM mDnaA ...
    mGcrA mSciP mCtrA CcrM DnaA GcrA SciP CtrA Sup DivKp I_DnaA I_GcrA I_ccrM tot;

global CPLX1 CpdR CpdRP CPLX2 RcdA CPLX3; 

global p;

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%&&&&&&&&%

% parameters for which cell we are tracking on
global CELLTYPE;

% cell type track
H = CELLTYPE; 
% parameters for the equations

% DNA replication, methylation
kelong = 0.00653;
nelong = 4;
Pelong = 0.05;
kmcori = 1.4;	Jmcori = 0.95;	nmcori = 4;
kmccrM = 1.4;	JmccrM = 0.95;	nmccrM = 4;
kmctrA = 1.4;	JmctrA = 0.95;	nmctrA = 4;

% DnaA, GcrA, CtrA, CcrM, SciP
ksDna = 0.1;   kdDna = 0.12;
ksGcrA = 0.045;  kdGcrA = 0.04;
ksCcrM = 0.1;  kdCcrM = 0.1;
ksSciP = 0.136;  kdSciP = 0.072;
ksCtrA = 0.036;  kdCtrA = 0.002;
kdCtrADivKp = 0.15; ndCtrADivKp = 2; JdCtrADivKp = 1;

ksI_DnaA = 0.05; kdI_DnaA = 0.05;
ksI_GcrA = 0.005; kdI_GcrA = 0.005;
ksI_ccrM = 0.1; kdI_ccrM = 0.1; % adjust by MC

%mCcrM
ksmCcrM = 0.256;     kdmCcrM = 0.08;  % adjust by MC
JaCcrMCtrA = 5;   naCcrMCtrA = 2;
JI_ccrMSciP = 6;   nI_ccrMSciP = 2;

%mDnaA
ksmDnaA = 0.0605;   kdmDnaA = 0.015; % adjust by MC
JI_DnaAGcrA = 3;	nI_DnaAGcrA = 2;
JaDnaACtrA = 5;   naDnaACtrA = 2;

%mGcrA
ksmGcrA = 5.6;    kdmGcrA = 0.6; % adjust by MC
JI_GcrACtrA = 5;   nI_GcrACtrA = 2;
JaGcrADnaA = 1.25;   naGcrADnaA = 2;

%mSciP
ksmSciP = 0.5;    kdmSciP = 0.06; % adjust by MC
JaSciPCtrA = 5; naSciPCtrA = 2;

%mCtrA
ks1mCtrA = 0.09; 
ks2mCtrA = 0.99;
kdmCtrA = 0.1; % adjust by MC
JaCtrACtrA = 5;	naCtrACtrA = 2;
JaCtrAGcrA = 3;	naCtrAGcrA = 2;
JiCtrACtrA = 8;	niCtrACtrA = 4;
JiCtrASciP = 8;	niCtrASciP = 4;

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
kaIni = 0.05;	kIIni = 1; JaIni = 1; naIni = 4; JIIni = 1; nIIni = 4; 
% kaIni = 0.05; % adjust by MC
thetaCtrA = 1;	nthetaCtrA = 4;
thetaDnaA = 3;	nthetaDnaA = 4;

dydt(Ini) = kIIni*(JIIni^nIIni/(JIIni^nIIni + (y(CtrA)/thetaCtrA)^nthetaCtrA) ...
     * kaIni*(y(DnaA)/thetaDnaA)^nthetaDnaA/(JaIni^naIni + (y(DnaA)/thetaDnaA)^nthetaDnaA) );

dydt(Elong) = kelong*y(Elong)^nelong/(y(Elong)^nelong + Pelong^nelong)*y(Count);

dydt(DNA) = kelong*y(Elong)^nelong/(y(Elong)^nelong + 0.05^nelong)*y(Count);

dydt(Count) = 0;

% DNA methylation
dydt(hcori) =  - kmcori*y(CcrM)^nmcori/(Jmcori^nmcori + y(CcrM)^nmcori)*y(hcori);
dydt(hCcrM) =  - kmccrM*y(CcrM)^nmctrA/(JmctrA^nmctrA + y(CcrM)^nmctrA)*y(hCcrM);
dydt(hCtrA) =  - kmctrA*y(CcrM)^nmccrM/(JmccrM^nmccrM + y(CcrM)^nmccrM)*y(hCtrA);

% mCcrM mDnaA mGcrA mSciP mCtrA
dydt(I_ccrM) = ksI_ccrM*(y(CtrA)^naCcrMCtrA/(JaCcrMCtrA^naCcrMCtrA + y(CtrA)^naCcrMCtrA) ...
              * JI_ccrMSciP^nI_ccrMSciP/(JI_ccrMSciP^nI_ccrMSciP + y(SciP)^nI_ccrMSciP) ) ...
              * y(hCcrM) - kdI_ccrM*y(I_ccrM);
          
dydt(mCcrM) = ksmCcrM*y(I_ccrM) - kdmCcrM*y(mCcrM);
      
dydt(mDnaA) = ( ksmDnaA*JI_DnaAGcrA^nI_DnaAGcrA/(JI_DnaAGcrA^nI_DnaAGcrA + y(GcrA)^nI_DnaAGcrA) ...
              * y(CtrA)^naDnaACtrA/(JaDnaACtrA^naDnaACtrA + y(CtrA)^naDnaACtrA) ) ...
              * (2 - y(hcori)) - kdmDnaA*y(mDnaA);
          
dydt(mGcrA) = ( ksmGcrA*y(DnaA)^naGcrADnaA/(JaGcrADnaA^naGcrADnaA + y(DnaA)^naGcrADnaA) ...
              * JI_GcrACtrA^nI_GcrACtrA/(JI_GcrACtrA^nI_GcrACtrA + y(CtrA)^nI_GcrACtrA) ) ...
              - kdmGcrA*y(mGcrA);

dydt(mSciP) = ksmSciP*y(CtrA)^naSciPCtrA/(JaSciPCtrA^naSciPCtrA + y(CtrA)^naSciPCtrA) ...
              - kdmSciP*y(mSciP);


dydt(mCtrA) = (ks1mCtrA*y(CtrA)^naCtrACtrA/(JaCtrACtrA^naCtrACtrA + y(CtrA)^naCtrACtrA)) ...
              + (ks2mCtrA*y(GcrA)^naCtrAGcrA/(JaCtrAGcrA^naCtrAGcrA + y(GcrA)^naCtrAGcrA) ...
              *JiCtrACtrA^niCtrACtrA/(JiCtrACtrA^niCtrACtrA + y(CtrA)^niCtrACtrA) ...
              *JiCtrASciP^niCtrASciP/(JiCtrASciP^niCtrASciP + y(SciP)^niCtrASciP) )...
              *y(hCtrA) - kdmCtrA*y(mCtrA) ;

% CcrM DnaA GcrA SciP CtrA CtrA
dydt(CcrM) = ksCcrM*y(mCcrM) - kdCcrM*y(CcrM);

dydt(I_DnaA) = ksI_DnaA*y(mDnaA) - kdI_DnaA*y(I_DnaA);

dydt(DnaA) = ksDna*y(mDnaA) - kdDna*y(DnaA); % adjust by MC

dydt(I_GcrA) = ksI_GcrA*y(mGcrA) - kdI_GcrA*y(I_GcrA);

dydt(GcrA) = ksGcrA*y(mGcrA) - kdGcrA*y(GcrA); % adjust by MC

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
%p.divkp=0.004218*t_d+ 0.5779;
%p.divkp(p.divkp<0)=0;
%0.92
% Linear model Poly3:
%  p.divkp = 1.64e-06*t_d.^3-0.0003057 *t_d.^2 + 0.02095*t_d +0.3282 ;
  
% Coefficients (with 95% confidence bounds):
%linear model poly1 divkp
p.divkp=0.004218*t_d+0.5779;
    
     
%RcdA
q1 =  -8.725e-08 ;
       q2 =   2.948e-05;
       q3 =   -0.003332;
       q4 =      0.1384 ;
       q5 =     -0.8569;
       
p.rcda=q1*t_d.^4 + q2*t_d.^3 + q3*t_d.^2 + q4*t_d+ q5;


%% cdG

cdG= 0.3233 *sin(pi/75*t_d+0.3469 )+ 0.3363;
cdG(cdG<0)=0;
%%
%p.divkp/(p.divkp+p.J2) this term was removed from CpdrP and CpdR
%8/30/2020.

dydt(CPLX1)=p.k1_pos*p.clpxp*y(CpdR)-p.k1_neg*y(CPLX1)...
    -p.k3_pos*y(CPLX1)*y(RcdA)+p.k3_neg*y(CPLX2);

dydt(CpdR)=p.ks_cpdr-p.kd_cpdr*y(CpdR)*y(CPLX1)/(y(CPLX1)+p.J1)...
    +p.k1_neg*y(CPLX1)-p.k1_pos*p.clpxp*y(CpdR)...
    +p.k2_pos*y(CpdRP)-p.k2_neg*y(CpdR);

dydt(CpdRP)=p.k2_neg*y(CpdR)-p.k2_pos*y(CpdRP);


dydt(CPLX2)=p.k3_pos*y(CPLX1)*y(RcdA)-p.k3_neg*y(CPLX2)+p.k5_neg*y(CPLX3)-p.k5_pos*cdG^2*y(CPLX2);
dydt(RcdA)=p.ks_rcda*y(RcdA)^1/(y(RcdA)^1+p.J3^1)-p.kd_rcda*y(RcdA)*y(CPLX1)/(y(CPLX1)+p.J4);


dydt(CPLX3)=-p.k5_neg*y(CPLX3)+p.k5_pos*cdG^2*y(CPLX2);

% end of equations
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%&&&&&&&&%

end
