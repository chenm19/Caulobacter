function dydt = ODEs5regCPLX(t, y)
%%variable map

global Ini Elong DNA Count hcori hCcrM hCtrA mCcrM mDnaA ...
    mGcrA mSciP mCtrA CcrM DnaA GcrA SciP CtrA Sup DivKp I II III tot;

global CPLX1 CpdR CpdRP CPLX2 RcdA CPLX3; 

global p;

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
ss = 1; sd = 1;
ksCcrM = ss*0.1;  kdCcrM = sd*0.1; % adjust by MC
ksDna = 0.1;   kdDna = 0.1;
ss = 1; sd = 1.1;
ksDna = ss*0.1;   kdDna = sd*0.1; % adjust by MC
ksGcrA = 0.1;  kdGcrA = 0.1;
ss = 0.5; sd =0.4;
ksGcrA = ss*0.1;  kdGcrA = sd*0.1;% adjust by MC
ksSciP = 0.08;  kdSciP = 0.08;
ss = 1.7; sd = .9;
ksSciP = ss*0.08;  kdSciP = sd*0.08; % adjust by MC
ksCtrA = 0.024;  kdCtrA = 0.002;
ss = 1.5; sd = 1;
ksCtrA = ss*0.024;  kdCtrA = sd*0.002; % adjust by MC
kdCtrADivKp = 0.15; ndCtrADivKp = 2; JdCtrADivKp = 1;

ksI = 0.05; kdI = 0.05;
ksII = 0.05; kdII = 0.05;
ss = 0.1; sd = 0.1;
ksII = ss*0.05; kdII = sd*0.05; % adjust by MC
ksIII = 0.1; kdIII = 0.1;
ss = 1; sd = 1;
ksIII = ss*0.1; kdIII = sd*0.1; % adjust by MC

%mCcrM
ksmCcrM = 0.32;     kdmCcrM = 0.08;
ss = .9; sd = 1;
ksmCcrM = ss*0.32;     kdmCcrM = sd*0.08;  % adjust by MC
JaCcrMCtrA = 5;   naCcrMCtrA = 2;
JiCcrMSciP = 6;   niCcrMSciP = 2;

%mDnaA
ksmDnaA = 0.055;   kdmDnaA = 0.015;
ss = 1.2; sd = 1;
ksmDnaA = ss*0.055;   kdmDnaA =sd*0.015; % adjust by MC
JiDnaAGcrA = 3;	niDnaAGcrA = 2;
JaDnaACtrA = 5;   naDnaACtrA = 2;

%mGcrA
ksmGcrA = 2.8;    kdmGcrA = 0.3;
ss = 2; sd = 2;
ksmGcrA = ss*2.8;    kdmGcrA = sd*0.3; % adjust by MC
JiGcrACtrA = 5;   niGcrACtrA = 2;
JaGcrADnaA = 1.5;   naGcrADnaA = 2;

%mSciP
ksmSciP = 0.5; kdmSciP = 0.06;
ss = 1; sd = 1;
ksmSciP = ss*0.5;    kdmSciP = sd*0.06; % adjust by MC
JaSciPCtrA = 5; naSciPCtrA = 2;

%mCtrA
ksmCtrA = 0.9;	kdmCtrA = 0.1;
ss = 1; sd = 1;
ksmCtrA = ss*0.9;	kdmCtrA = sd*0.1; % adjust by MC
JaCtrACtrA = 5;	naCtrACtrA = 2;
JaCtrAGcrA = 3;	naCtrAGcrA = 2;
JiCtrACtrA = 5;	niCtrACtrA = 2;
JiCtrASciP = 6;	niCtrASciP = 2;

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
kaIni = 0.05;	kiIni = 1; JaIni = 1; naIni = 4; JiIni = 1; niIni = 4; 
% kaIni = 0.05; % adjust by MC
thetaCtrA = 1;	nthetaCtrA = 4;
thetaDnaA = 3;	nthetaDnaA = 4;

dydt(Ini) = kiIni*(JiIni^niIni/(JiIni^niIni + (y(CtrA)/thetaCtrA)^nthetaCtrA) ...
     * kaIni*(y(DnaA)/thetaDnaA)^nthetaDnaA/(JaIni^naIni + (y(DnaA)/thetaDnaA)^nthetaDnaA) );

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

          ss = 0.3; % adjust by MC
dydt(mCtrA) = (ksmCtrA*y(CtrA)^naCtrACtrA/(JaCtrACtrA^naCtrACtrA + y(CtrA)^naCtrACtrA) ...
              + ksmCtrA*y(GcrA)^naCtrAGcrA/(JaCtrAGcrA^naCtrAGcrA + y(GcrA)^naCtrAGcrA) ...
              *JiCtrACtrA^niCtrACtrA/(JiCtrACtrA^niCtrACtrA + y(CtrA)^niCtrACtrA) ...
              *JiCtrASciP^niCtrASciP/(JiCtrASciP^niCtrASciP + y(SciP)^niCtrASciP) )...
              *y(hCtrA) - kdmCtrA*y(mCtrA) ;

% CcrM DnaA GcrA SciP CtrA CtrA
dydt(CcrM) = ksCcrM*y(mCcrM) - kdCcrM*y(CcrM);
dydt(I) = ksI*y(mDnaA) - kdI*y(I);
% dydt(DnaA) = ksDna*y(I) - kdDna*y(DnaA);
dydt(DnaA) = ksDna*y(mDnaA) - kdDna*y(DnaA); % adjust by MC

dydt(II) = ksII*y(mGcrA) - kdII*y(II);
dydt(GcrA) = ksGcrA*y(II) - kdGcrA*y(GcrA);
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
