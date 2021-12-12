%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% Script for ordinary differential equations of Caulobacter cell cycle model
% Change the parameters to generate WILD type and all MUTANT cases.
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%

function dydt = odes5RegCPLX(t, y, para)

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

p.kmcori =  para(1);
p.kmccrM =  para(2);	
p.kmctrA =  para(3);
p.ksI_ccrM =  para(4);
p.kdI_ccrM =  para(5);
p.ksmCcrM =  para(6);
p.kdmCcrM =  para(7);
p.ksmDnaA =  para(8);
p.kdmDnaA =  para(9);
p.ksmGcrA = para(10);
p.kdmGcrA =  para(11);
p.ksmSciP =  para(12);
p.kdmSciP =  para(13);
p.ks1mCtrA =  para(14);
p.ks2mCtrA =  para(15);
p.kdmCtrA =  para(16);
p.ksDnaA =  para(17);
p.ksGcrA =  para(18);
p.ksCcrM =  para(19);
p.ksSciP =  para(20);
p.kdSciP =  para(21);
p.ksCtrA =  para(22);
p.kdCtrAClpXP =  para(23);
p.k1_pos=  para(24);
p.k1_neg=  para(25);
p.ks_cpdr=  para(26);
p.kd_cpdr=  para(27);
p.dephos_CpdR=  para(28);
p.phos_CpdR=  para(29);
p.k2_pos=    para(30);
p.k2_neg=    para(31);
p.ks_rcda=   para(32);
p.kd_rcda=   para(33);
p.k3_pos=    para(34);
p.k3_neg=    para(35);
p.kscdG=     para(36);
p.kdcdG=     para(37);
p.ksPleD=    para(38);
p.kdPleD=    para(39);
p.phosPleD=   para(40);   
p.dephosPleD=  para(41);
p.ksPdeA=   para(42);
p.kdPdeA=  para(43);
p.kphoCtrA = para(44);
p.kdephoCckA = para(45);
p.kdephoCtrA = para(46);
p.kphoCckA = para(47);


K3 = p.k3_neg/p.k3_pos;
fun = @(x) 4*x^3 - 4*(y(CPLX2)+y(cdG))*x^2 + (y(cdG)^2+4*y(cdG)*y(CPLX2)+K3)*x - y(cdG)^2*y(CPLX2);
options = optimoptions(@lsqnonlin,'display','off');
xsol = lsqnonlin(fun, 2, 0, min(y(cdG),y(CPLX2)), options);
y(CPLX3) = xsol;
%% mutant list%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p.ksmDnaA = 0; %%%%% DEL DnaA Mutant %%%%%%%%
% p.ksmCcrM = 0; %%%%%%%%% DEL CcrM Mutant %%%%%%%%%%%%%
% p.ksmGcrA = 0; %%%%%%%%% DEL GcrA Mutant %%%%%%%%%%%%%
% p.ks1mCtrA = 0;  %%%SM921 Mutant %%%%%%%%
% p.kdCtrAClpXP = 0.1*p.kdCtrAClpXP; %%%%%%% CtrA3 Mutant %%%%%%%
% p.kscdG = 0; %%%%%%%%%%% cdG MUTANT %%%%%%%%%%%%%%%%
% p.ksPleD = 0; %%%%%% PleD Mutant %%%%%%
% p.ksPdeA = 0; %%%%%% PdeA Mutant %%%%%%

% p.JdCtrAClpXP = 0;% prediction 1. CtrA deg mutant: no cyclic degradation of CtrA
% p.JdCpdR=0;%prediction 2. CpdR deg mutant: no cyclic degradation of CpdR
% p.JdRcdA=0;%prediction 3. RcdA deg mutant: no cyclic degradation of RcdA
% p.JdCtrAClpXP = 0; p.JdCpdR=0; p.JdRcdA=0;%prediction 4. no cyclic degradation in the system
dephoConstant=0; p5=1;
% p.JdCtrAClpXP=0; dephoConstant=5; p5=0;%prediction 5. deleting cdG impact on both cyclic degradation of CtrA and dephosophorylation of CckA

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%&&&&&&&&%
% Differential equations of the model
dydt = zeros(tot, 1);  

% Initiation
dydt(Ini) = p.ksIni * (1+1/(p.JiIni^4+(y(hcori)/p.thetaCori).^4)) ...
           *( (y(DnaA)/p.thetaDnaA)^4 / (p.JaIni^4 + (y(CtrAP)/p.thetaCtrA)^4 + (y(DnaA)/p.thetaDnaA)^4) ); 


% DNA replication
dydt(Elong) = p.kelong*y(Elong)^p.nelong/(y(Elong)^p.nelong + p.Pelong^p.nelong)*y(Count);
dydt(DNA) = p.kelong*y(Elong)^p.nelong/(y(Elong)^p.nelong + p.Pelong^p.nelong)*y(Count);
dydt(Count) = 0;


% DNA methylation 
dydt(hcori) =  - p.kmcori*y(CcrM)^p.nmcori/(p.Jmcori^p.nmcori + y(CcrM)^p.nmcori)*y(hcori);
dydt(hCcrM) =  - p.kmccrM*y(CcrM)^p.nmctrA/(p.JmctrA^p.nmctrA + y(CcrM)^p.nmctrA)*y(hCcrM);
dydt(hCtrA) =  - p.kmctrA*y(CcrM)^p.nmccrM/(p.JmccrM^p.nmccrM + y(CcrM)^p.nmccrM)*y(hCtrA);


%mCcrM

dydt(I_ccrM) = p.ksI_ccrM*(y(CtrAP)^p.naCcrMCtrA/(p.JaCcrMCtrA^p.naCcrMCtrA + y(CtrAP)^p.naCcrMCtrA) ...
              * p.JI_ccrMSciP^p.nI_ccrMSciP/(p.JI_ccrMSciP^p.nI_ccrMSciP + y(SciP)^p.nI_ccrMSciP) ) ...
              * y(hCcrM) - p.kdI_ccrM*y(I_ccrM);
dydt(mCcrM) = p.ksmCcrM*y(I_ccrM) - p.kdmCcrM*y(mCcrM);
      

%mDnaA
 
% dydt(mDnaA) = ( ksmDnaA*JI_DnaAGcrA^nI_DnaAGcrA/(JI_DnaAGcrA^nI_DnaAGcrA + y(GcrA)^nI_DnaAGcrA) ...
%               * y(CtrAP)^naDnaACtrA/(JaDnaACtrA^naDnaACtrA + y(CtrAP)^naDnaACtrA) ) ...
%               * (2 - y(hcori)) - kdmDnaA*y(mDnaA); %old mDnaA equation
          
dydt(mDnaA) =  p.ksmDnaA*p.JI_DnaAGcrA^p.nI_DnaAGcrA/(p.JI_DnaAGcrA^p.nI_DnaAGcrA + y(GcrA)^p.nI_DnaAGcrA) ...
              * (2 - y(hcori)) - p.kdmDnaA*y(mDnaA);
               
          
%mGcrA
% ksmGcrA = 5.6;    kdmGcrA = 0.6; 
% JI_GcrACtrA = 5;   nI_GcrACtrA = 2;
% JaGcrADnaA = 1.25;   naGcrADnaA = 2;

 
dydt(mGcrA) = ( p.ksmGcrA*y(DnaA)^p.naGcrADnaA/(p.JaGcrADnaA^p.naGcrADnaA + y(DnaA)^p.naGcrADnaA) ...
              * p.JI_GcrACtrA^p.nI_GcrACtrA/(p.JI_GcrACtrA^p.nI_GcrACtrA + y(CtrAP)^p.nI_GcrACtrA) ) ...
              - p.kdmGcrA*y(mGcrA);

%mSciP
% ksmSciP = 0.5;    kdmSciP = 0.06/1.5; 
% JaSciPCtrA = 5; naSciPCtrA = 2;

dydt(mSciP) = p.ksmSciP*y(CtrAP)^p.naSciPCtrA/(p.JaSciPCtrA^p.naSciPCtrA + y(CtrAP)^p.naSciPCtrA) ...
              - p.kdmSciP*y(mSciP);


%mCtrA

dydt(mCtrA) = (p.ks2mCtrA*y(CtrAP)^p.naCtrACtrA/(p.JaCtrACtrA^p.naCtrACtrA + y(CtrAP)^p.naCtrACtrA)) ...
              + (p.ks1mCtrA*y(GcrA)^p.naCtrAGcrA/(p.JaCtrAGcrA^p.naCtrAGcrA + y(GcrA)^p.naCtrAGcrA) ...
              *p.JiCtrACtrA^p.niCtrACtrA/(p.JiCtrACtrA^p.niCtrACtrA + y(CtrAP)^p.niCtrACtrA) ...
              *p.JiCtrASciP^p.niCtrASciP/(p.JiCtrASciP^p.niCtrASciP + y(SciP)^p.niCtrASciP) )...
              *y(hCtrA) - p.kdmCtrA*y(mCtrA) ;

           
  

%%new 1/26/2021 RX 
 
dydt(CcrM) = p.ksCcrM*y(mCcrM) - p.kdCcrM*y(CcrM);
% dydt(I_DnaA) = ksI_DnaA*y(mDnaA) - kdI_DnaA*y(I_DnaA);
dydt(DnaA) = p.ksDnaA*y(mDnaA) - p.kdDnaA*y(DnaA); % adjust by MC
% dydt(I_GcrA) = ksI_GcrA*y(mGcrA) - kdI_GcrA*y(I_GcrA);
dydt(GcrA) = p.ksGcrA*y(mGcrA) - p.kdGcrA*y(GcrA); % adjust by MC
dydt(SciP) = p.ksSciP*y(mSciP) - p.kdSciP*y(SciP);
dydt(CtrA) = p.ksCtrA*y(mCtrA) - (p.kdCtrA + p.kdCtrAClpXP*y(CPLX3)^p.ndCtrAClpXP/...
             (p.JdCtrAClpXP^p.ndCtrAClpXP + y(CPLX3)^p.ndCtrAClpXP) )*y(CtrA) ...
             -p.kphoCtrA*y(CckAP)*y(CtrA)+p.kdephoCtrA*y(CtrAP);%new 1/26/2021 RX

% Zring
dydt(Zring) = ksZring;

% DivKp
dydt(DivKp) = 0;

%%


T=150;        % period of Caulobacter
t_d=rem(t,T); % return remainder after division t/T

  
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
    -p.k2_pos*y(CPLX1)*y(RcdA)+p.k2_neg*(y(CPLX2)-y(CPLX3));

dydt(CpdR)=p.ks_cpdr*y(CtrAP)^p.naCpdRCtrA/(p.JaCpdRCtrA^p.naCpdRCtrA+y(CtrAP)^p.naCpdRCtrA)-p.kd_cpdr*y(CpdR)*y(CPLX1)/(y(CPLX1)+p.JdCpdR)...
    +p.k1_neg*y(CPLX1)-p.k1_pos*p.clpxp*y(CpdR)...
    +p.dephos_CpdR*y(CpdRP)-p.phos_CpdR*y(CpdR)*y(CckAP);%1/26/2021 RX
% dydt(CpdR)=p.ks_cpdr-p.kd_cpdr*y(CpdR)*y(CPLX1)/(y(CPLX1)+p.J1)...
%     +p.k1_neg*y(CPLX1)-p.k1_pos*p.clpxp*y(CpdR)...
%     +p.k2_pos*y(CpdRP)-p.k2_neg*y(CpdR)*y(CckAP);%2/6/2021 RX

dydt(CpdRP)=-p.kd_cpdr*y(CpdRP)*y(CPLX1)/(y(CPLX1)+p.JdCpdR)+p.phos_CpdR*y(CpdR)*y(CckAP)-p.dephos_CpdR*y(CpdRP);%1/26/2021 RX

dydt(CPLX2)=p.k2_pos*y(CPLX1)*y(RcdA)-p.k2_neg*(y(CPLX2)-y(CPLX3));%+p.k3_neg*y(CPLX3)-p.k3_pos*y(cdG)^2*y(CPLX2);
dydt(RcdA)=0+p.ks_rcda*y(CtrAP)^p.naRcdACtrA/(y(CtrAP)^p.naRcdACtrA+p.JaRcdACtrA^p.naRcdACtrA)-p.kd_rcda*y(RcdA)*y(CPLX1)/(y(CPLX1)+p.JdRcdA);
%1/26/2021 RX
dydt(CPLX3)=0;

%%new % 1/26/2021 RX
dydt(CtrAP) = - (p.kdCtrA + p.kdCtrAClpXP*y(CPLX3)^p.ndCtrAClpXP/...
             (p.JdCtrAClpXP^p.ndCtrAClpXP + y(CPLX3)^p.ndCtrAClpXP) )*y(CtrAP)...
             +p.kphoCtrA*y(CckAP)*y(CtrA)-p.kdephoCtrA*y(CtrAP);

dydt(CckAP) = p.kphoCckA*(p.CckAT-y(CckAP))-p.kdephoCckA*(1+p5*p.alpha_cdG*(y(cdG)-2*y(CPLX3))+dephoConstant)*y(CckAP);

%% cdG
% dydt(cdG)=p.kscdG*(1+p.alpha_PleD*y(PleDP))*p.JicdGcdG^2/(p.JicdGcdG^2+y(cdG)^2)...
%     -p.kdcdG*(1+p.alpha_PdeA*y(PdeA))*y(cdG); %/(y(cdG)+p.JdcdG);

dydt(cdG)=p.kscdG*(1+p.alpha_PleD*y(PleDP))*p.JicdGcdG^2/(p.JicdGcdG^2+(y(cdG)-2*y(CPLX3))^2)...
    -p.kdcdG*(1+p.alpha_PdeA*y(PdeA))*(y(cdG)-2*y(CPLX3));
dydt(PleD)=p.ksPleD*y(CtrAP)^p.naPleDCtrA/(y(CtrAP)^p.naPleDCtrA+p.JaPleDCtrA^p.naPleDCtrA)...
    -p.kdPleD*y(PleD)-p.phosPleD*y(PleD)+p.dephosPleD*p.PleC*y(PleDP);

dydt(PdeA)=p.ksPdeA*y(CtrAP)^p.naPdeACtrA/(y(CtrAP)^p.naPdeACtrA+p.JaPdeACtrA^p.naPdeACtrA)...
    -p.kdPdeA*y(PdeA)*y(CPLX1)/(y(CPLX1)+p.JdPdeA);
dydt(PleDP)=p.phosPleD*y(PleD)-p.dephosPleD*p.PleC*y(PleDP);

% end of equations
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%&&&&&&&&%

end