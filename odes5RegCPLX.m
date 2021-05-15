function dydt = odes5RegCPLX(t, y)
%%variables
global Ini Elong DNA Count hcori hCcrM hCtrA ...
       mCcrM mDnaA mGcrA mSciP mCtrA ...
       CcrM DnaA GcrA SciP CtrA Zring DivKp ...
       I_DnaA I_GcrA I_ccrM tot;

global CPLX1 CpdR CpdRP CPLX2 RcdA CPLX3; 

global CtrAP CckAP cdG PleD PdeA; %new 1/26/2021 RX
global PleDP; %13/3/2021 RX
global p;
global ksZring;

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%&&&&&&&&%
% parameters for which cell we are tracking on
global CELLTYPE;
% cell type track
H = CELLTYPE; 

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%&&&&&&&&%
% Differential equations of the model
dydt = zeros(tot, 1);  

% Initiation
ksIni = 0.00031044; JiIni = 1.4565; JaIni = 1;
thetaDnaA = 0.5; thetaCtrA = 6; thetaCori = 0.308;

dydt(Ini) = ksIni * (1+1/(JiIni^4+(y(hcori)/thetaCori).^4)) ...
            *( (y(DnaA)/thetaDnaA)^4 / (JaIni^4 + (y(CtrAP)/thetaCtrA)^4 + (y(DnaA)/thetaDnaA)^4) ); 

% dydt(Ini) = 0.00071*(y(DnaA)/0.5)^4 / (1 + (y(CtrAP)/6)^4 + (y(DnaA)/0.5)^4) * ...
%      (1.1-y(hcori).^4./(0.05+y(hcori).^4)/5)/2.06; 
 
% dydt(Ini) = 0.001; % fake init to tune parameters

% kIIni = 1; kaIni = 0.05*20; 
% JaIni = 1; naIni = 4; JIIni = 1; nIIni = 4; 
% thetaCtrA = 0.1;	nthetaCtrA = 2;
% thetaDnaA = 1;	nthetaDnaA = 2;

% dydt(Ini) = kIIni*(JIIni^nIIni/(JIIni^nIIni + (y(CtrAP)/thetaCtrA)^nthetaCtrA) ...
%      * kaIni*(y(DnaA)/thetaDnaA)^nthetaDnaA/(JaIni^naIni + (y(DnaA)/thetaDnaA)^nthetaDnaA) );
 
% nthetaCtrA = 4; nthetaDnaA = 4; 
% kaIni = 0.0022;  %2/18/21 MC
% JaIni = 1;
% thetaDnaA = 2; thetaCtrA = 6; thetaCori = 0.0001;
% dydt(Ini) = kaIni*(y(DnaA)/thetaDnaA)^nthetaDnaA / ...
%      (JaIni + (y(CtrAP)/thetaCtrA)^nthetaCtrA  ...
%         + (y(DnaA)/thetaDnaA)^nthetaDnaA ...
%         + (y(hcori)/thetaCori));
  
% DNA replication
kelong = 0.00653; nelong = 4; Pelong = 0.05;

dydt(Elong) = kelong*y(Elong)^nelong/(y(Elong)^nelong + Pelong^nelong)*y(Count);
dydt(DNA) = kelong*y(Elong)^nelong/(y(Elong)^nelong + Pelong^nelong)*y(Count);
dydt(Count) = 0;


% DNA methylation 
kmcori = 1.4;	Jmcori = 0.95;	nmcori = 4; 
kmccrM = 1.4;	JmccrM = 0.95;	nmccrM = 4;
kmctrA = 1.4;	JmctrA = 0.95;	nmctrA = 4;

dydt(hcori) =  - kmcori*y(CcrM)^nmcori/(Jmcori^nmcori + y(CcrM)^nmcori)*y(hcori);
dydt(hCcrM) =  - kmccrM*y(CcrM)^nmctrA/(JmctrA^nmctrA + y(CcrM)^nmctrA)*y(hCcrM);
dydt(hCtrA) =  - kmctrA*y(CcrM)^nmccrM/(JmccrM^nmccrM + y(CcrM)^nmccrM)*y(hCtrA);


%mCcrM
ksI_ccrM = 0.1; kdI_ccrM = 0.1/1.5; 
ksmCcrM = 0.256;  kdmCcrM = 0.08; 
JaCcrMCtrA = 5;   naCcrMCtrA = 2;
JI_ccrMSciP = 6;   nI_ccrMSciP = 2;

%ksmCcrM = 0; %%%%%%%%% DEL CcrM Mutant %%%%%%%%%%%%%

dydt(I_ccrM) = ksI_ccrM*(y(CtrAP)^naCcrMCtrA/(JaCcrMCtrA^naCcrMCtrA + y(CtrAP)^naCcrMCtrA) ...
              * JI_ccrMSciP^nI_ccrMSciP/(JI_ccrMSciP^nI_ccrMSciP + y(SciP)^nI_ccrMSciP) ) ...
              * y(hCcrM) - kdI_ccrM*y(I_ccrM);
dydt(mCcrM) = ksmCcrM*y(I_ccrM) - kdmCcrM*y(mCcrM);
      

%mDnaA
ksmDnaA = 0.0605*4;   kdmDnaA = 0.015*4; 
JI_DnaAGcrA = 3;	nI_DnaAGcrA = 2;
% JaDnaACtrA = 5;   naDnaACtrA = 2;

% ksmDnaA = 0; %%%%% DEL DnaA Mutant %%%%%%%%
 
% dydt(mDnaA) = ( ksmDnaA*JI_DnaAGcrA^nI_DnaAGcrA/(JI_DnaAGcrA^nI_DnaAGcrA + y(GcrA)^nI_DnaAGcrA) ...
%               * y(CtrAP)^naDnaACtrA/(JaDnaACtrA^naDnaACtrA + y(CtrAP)^naDnaACtrA) ) ...
%               * (2 - y(hcori)) - kdmDnaA*y(mDnaA); %old mDnaA equation
          
dydt(mDnaA) =  ksmDnaA*JI_DnaAGcrA^nI_DnaAGcrA/(JI_DnaAGcrA^nI_DnaAGcrA + y(GcrA)^nI_DnaAGcrA) ... ...
              * (2 - y(hcori)) - kdmDnaA*y(mDnaA);
               
          
%mGcrA
ksmGcrA = 5.6;    kdmGcrA = 0.6; 
JI_GcrACtrA = 5;   nI_GcrACtrA = 2;
JaGcrADnaA = 1.25;   naGcrADnaA = 2;

% ksmGcrA = 0; %%%%%%%%% DEL GcrA Mutant %%%%%%%%%%%%%
 
dydt(mGcrA) = ( ksmGcrA*y(DnaA)^naGcrADnaA/(JaGcrADnaA^naGcrADnaA + y(DnaA)^naGcrADnaA) ...
              * JI_GcrACtrA^nI_GcrACtrA/(JI_GcrACtrA^nI_GcrACtrA + y(CtrAP)^nI_GcrACtrA) ) ...
              - kdmGcrA*y(mGcrA);

%mSciP
ksmSciP = 0.5;    kdmSciP = 0.06/1.5; 
JaSciPCtrA = 5; naSciPCtrA = 2;

dydt(mSciP) = ksmSciP*y(CtrAP)^naSciPCtrA/(JaSciPCtrA^naSciPCtrA + y(CtrAP)^naSciPCtrA) ...
              - kdmSciP*y(mSciP);


%mCtrA
ks1mCtrA = 0.99; ks2mCtrA = 0.09;  kdmCtrA = 0.1/1.2; 
JaCtrACtrA = 5;	naCtrACtrA = 2;
JaCtrAGcrA = 3;	naCtrAGcrA = 2;
JiCtrACtrA = 8;	niCtrACtrA = 4;
JiCtrASciP = 8;	niCtrASciP = 4;

% ks1mCtrA = 0;  %%%SM921 Mutant %%%%%%%%

dydt(mCtrA) = (ks2mCtrA*y(CtrAP)^naCtrACtrA/(JaCtrACtrA^naCtrACtrA + y(CtrAP)^naCtrACtrA)) ...
              + (ks1mCtrA*y(GcrA)^naCtrAGcrA/(JaCtrAGcrA^naCtrAGcrA + y(GcrA)^naCtrAGcrA) ...
              *JiCtrACtrA^niCtrACtrA/(JiCtrACtrA^niCtrACtrA + y(CtrAP)^niCtrACtrA) ...
              *JiCtrASciP^niCtrASciP/(JiCtrASciP^niCtrASciP + y(SciP)^niCtrASciP) )...
              *y(hCtrA) - kdmCtrA*y(mCtrA) ;

          
% CcrM DnaA GcrA SciP CtrA CtrA
ksSciP = 0.136/1.15;  kdSciP = 0.072/1.2; %%2/18/2021 MC
ksDnaA = 0.0065*10; kdDnaA = 0.007*10; % degradation rates from literature, adjusted, otherwise too slow
ksGcrA = 0.028;  kdGcrA = 0.022; % degradation rates from literature
ksCcrM = 0.085;  kdCcrM = 0.07;  % degradation rates from literature
ksCtrA = 0.108/2.5; kdCtrA = 0.002; kdCtrAClpXP = 0.15*0.4; % degradation rates from literature, kdCtrADivKp in range of 0.5-1.7
ndCtrAClpXP= 2; JdCtrAClpXP = 1*4;
 
% kdCtrADivKp = 0.1*kdCtrADivKp; %%%%%%% CtrA3 Mutant %%%%%%%
  

%%new 1/26/2021 RX
CckAT=0.3;
kdephoCtrA =0.1;% log(2)/5;
kphoCtrA = 5;%7*kdephoCtrA/(93*0.9*CckAT);
kphoCckA = 1; kdephoCckA = 1; alpha = 10; 
 
dydt(CcrM) = ksCcrM*y(mCcrM) - kdCcrM*y(CcrM);
% dydt(I_DnaA) = ksI_DnaA*y(mDnaA) - kdI_DnaA*y(I_DnaA);
dydt(DnaA) = ksDnaA*y(mDnaA) - kdDnaA*y(DnaA); % adjust by MC
% dydt(I_GcrA) = ksI_GcrA*y(mGcrA) - kdI_GcrA*y(I_GcrA);
dydt(GcrA) = ksGcrA*y(mGcrA) - kdGcrA*y(GcrA); % adjust by MC
dydt(SciP) = ksSciP*y(mSciP) - kdSciP*y(SciP);
dydt(CtrA) = ksCtrA*y(mCtrA) - (kdCtrA + kdCtrAClpXP*y(CPLX3)^ndCtrAClpXP/...
             (JdCtrAClpXP^ndCtrAClpXP + y(CPLX3)^ndCtrAClpXP) )*y(CtrA) ...
             -kphoCtrA*y(CckAP)*y(CtrA)+kdephoCtrA*y(CtrAP);%new 1/26/2021 RX

% Zring
dydt(Zring) = ksZring;

% DivKp
dydt(DivKp) = 0;

%%
% p.kscdG = 0; %%%%%%%%%%% cdG MUTANT %%%%%%%%%%%%%%%%
% p.ksPleD = 0; %%%%%% PleD Mutant %%%%%%
% p.ksPdeA = 0; %%%%%% PdeA Mutant %%%%%%


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
    

%p.divkp/(p.divkp+p.J2) this term was removed from CpdrP and CpdR
%8/30/2020.
%PleC 3/13/2021 RX
       a1 =       80.09;
       b1 =     0.01299;
       c1 =        1.74;
       a2 =       78.77;
       b2 =     0.01333;
       c2 =       4.854;
p.PleC =  a1*sin(b1*t_d+c1) + a2*sin(b2*t_d+c2);
%%new 1/26/2021 RX
dydt(CPLX1)=p.k1_pos*p.clpxp*y(CpdR)-p.k1_neg*y(CPLX1)...
    -p.k2_pos*y(CPLX1)*y(RcdA)+p.k2_neg*y(CPLX2);

dydt(CpdR)=p.ks_cpdr*y(CtrAP)^p.naCpdRCtrA/(p.JaCpdRCtrA^p.naCpdRCtrA+y(CtrAP)^p.naCpdRCtrA)-p.kd_cpdr*y(CpdR)*y(CPLX1)/(y(CPLX1)+p.JdCpdR)...
    +p.k1_neg*y(CPLX1)-p.k1_pos*p.clpxp*y(CpdR)...
    +p.dephos_CpdR*y(CpdRP)-p.phos_CpdR*y(CpdR)*y(CckAP);%1/26/2021 RX
% dydt(CpdR)=p.ks_cpdr-p.kd_cpdr*y(CpdR)*y(CPLX1)/(y(CPLX1)+p.J1)...
%     +p.k1_neg*y(CPLX1)-p.k1_pos*p.clpxp*y(CpdR)...
%     +p.k2_pos*y(CpdRP)-p.k2_neg*y(CpdR)*y(CckAP);%2/6/2021 RX

dydt(CpdRP)=-p.kd_cpdr*y(CpdRP)*y(CPLX1)/(y(CPLX1)+p.JdCpdR)+p.phos_CpdR*y(CpdR)*y(CckAP)-p.dephos_CpdR*y(CpdRP);%1/26/2021 RX

dydt(CPLX2)=p.k2_pos*y(CPLX1)*y(RcdA)-p.k2_neg*y(CPLX2)+p.k3_neg*y(CPLX3)-p.k3_pos*y(cdG)^2*y(CPLX2);
dydt(RcdA)=0+p.ks_rcda*y(CtrAP)^p.naRcdACtrA/(y(CtrAP)^p.naRcdACtrA+p.JaRcdACtrA^p.naRcdACtrA)-p.kd_rcda*y(RcdA)*y(CPLX1)/(y(CPLX1)+p.JdRcdA);
% dydt(RcdA)=p.ks_rcda-p.kd_rcda*y(RcdA)*y(CPLX1)/(y(CPLX1)+p.J4);
%1/26/2021 RX
dydt(CPLX3)=-p.k3_neg*y(CPLX3)+p.k3_pos*y(cdG)^2*y(CPLX2);

%%new % 1/26/2021 RX
dydt(CtrAP) = - (kdCtrA + kdCtrAClpXP*y(CPLX3)^ndCtrAClpXP/...
             (JdCtrAClpXP^ndCtrAClpXP + y(CPLX3)^ndCtrAClpXP) )*y(CtrAP)...
             +kphoCtrA*y(CckAP)*y(CtrA)-kdephoCtrA*y(CtrAP);

dydt(CckAP) = kphoCckA*(CckAT-y(CckAP))-kdephoCckA*(1+alpha*y(cdG))*y(CckAP);
%% cdG
dydt(cdG)=p.kscdG*(1+p.DGC*y(PleDP))*p.JicdGcdG^2/(p.JicdGcdG^2+y(cdG)^2)...
    -p.kdcdG*(1+p.PDE*y(PdeA))*y(cdG); %/(y(cdG)+p.JdcdG);
% dydt(cdG)=p.kscdG*(1+p.DGC*y(PleD)*p.JicdGcdG^2/(p.JicdGcdG^2+y(cdG)^2))...
%     -p.kdcdG*(1+p.PDE*y(PdeA))*y(cdG)/(y(cdG)+p.JdcdG);

dydt(PleD)=p.ksPleD*y(CtrAP)^p.naPleDCtrA/(y(CtrAP)^p.naPleDCtrA+p.JaPleDCtrA^p.naPleDCtrA)...
    -p.kdPleD*y(PleD)-p.phosPleD*y(PleD)+p.dephosPleD*p.PleC*y(PleDP);

dydt(PdeA)=p.ksPdeA*y(CtrAP)^p.naPdeACtrA/(y(CtrAP)^p.naPdeACtrA+p.JaPdeACtrA^p.naPdeACtrA)...
    -p.kdPdeA*y(PdeA)*y(CPLX1)/(y(CPLX1)+p.JdPdeA);
dydt(PleDP)=p.phosPleD*y(PleD)-p.dephosPleD*p.PleC*y(PleDP);

% end of equations
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%&&&&&&&&%

end
