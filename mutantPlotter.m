%graph cellcycle
function [] = FinalFigGenerator(tout, yout)

set(groot, 'DefaultAxesFontSize', 28)
%set(groot, 'DefaultFigurePosition', [100 100 500 400])
set(groot,'DefaultFigureColormap',jet)  
% set(groot,'DefaultAxesColorOrder',[0 0 1; 0 0 1; ])
set(groot,'DefaultFigureGraphicsSmoothing','off')

Ini = 1;
Elong = 2;
DNA = 3;
Count = 4;
hcori = 5;
hCcrM = 6;
hCtrA = 7;
mCcrM = 8;
mDnaA = 9;
mGcrA = 10;
mSciP = 11;
mCtrA = 12;
CcrM = 13;
DnaA = 14;
GcrA = 15;
SciP = 16;
CtrA = 17;
Sup = 18;
DivKp = 19;
I = 20;
II = 21;
III = 22;
tot = 33;

CPLX1 = 23;
CpdR = 24;
CpdRP = 25;
CPLX2 = 26;
RcdA = 27;
CPLX3 = 28;

CtrAP = 29;
CckAP=30; 
cdG=31; 
PleD=32; 
PdeA=33;
%data
t =  [5, 30, 60, 90, 120, 150];
dDnaA = [1262	1164 500	496 	1066	1028]./445;
dCcrM = [17	17	24	445	435	101]./445;
dSciP = [2858	451	199	796	2956	2649]./445;
dGcrA = [550	3001	1275	476	725	1493]./445;
dCtrA = [195	612	2251	3216	2486	1617]./445;
% publication 1:
tp1 = [0 20 40 60 80 100 120 140]  ;
pDnaA1 = [4784.305 10537.962 9511.012 3134.012 1447.184 1014.355 3693.083 6339.548];
pGcrA1 = [0 8632.062 16136.355 14185.012 14835.891 8353.184 9923.598 10884.719];
pCtrA1 = [6648.012 0 0 0 4042.598 8955.376 7497.598 6225.962];

%publication 2:
tp2 = [0 20 40 60 80 100 120 140]  ;
pSciP2 = [9239.255 9518.134 4585.305 610.042 161.657 388.435 4424.062 10748.426];
pGcrA2 = [0 444.87 4548.648 6148.669 5947.891 4592.062 5067.406 5809.548];
pCtrA2 = [7177.962 9270.255 1045.406 0 315.506 5709.77 8172.598 7724.355];

%publication 3:
tp3 = [0 16.5 33 49.5 66 82.5 99 115.5 132 148.5]  ;
pCcrM3 = [1314.983 1261.619 2128.518 1517 2385.477 2660.548 2990.276 7318.69 7300.154 4928.518];
pCtrA3 = [4366.669 544.698 678.284 1565.962 3814.74 6857.891 7663.376 13218.74 12383.598 10310.083];

%publication 4:
tp4 = [0 10 20 30 40 50 60 70 80 90 100]  ;
pSciP4 = [8176.062 3082.648 894.456 248.12 248.12 248.12 248.12 855.77 1727.527 2709.941 2193.042];
pCtrA4 = [3465.406 1278.77 928.527 590.991 1960.355 3866.477 4892.77 5509.062 4213.941 4028.77 2851.042];

%publication 5:
tp5 = [0 20 40 60 80 100 120 140]  ;
pDnaA5 = [6925.376 14430.033 13141.619 4829.376 2257.719 2059.497 6289.912 9555.912];
pGcrA5 = [0 8417.062 16478.891 14706.962 16155.255 9149.255 10833.548 12333.669];
pCtrA5 = [8860.134 0 0 0 4023.77 12430.841 10419.669 7582.255];
pCcrM5 = [0 0 0 0 166 3633.548 10993.962 13607.912];

pCpdRP = [9118.983,4425.406,3905.163,1912.92,5768.719,8952.497,6521.347];
tpCpdRP = [20,40,60,80,100,120,140]  ;

pCpdR = [17942.489,19085.196,12028.882,9405.711,9145.711,9438.782,12948.823];
tpCpdR = [0, 20,40,60,80,100,120]  ;

rcda =[0.0000136,0.093,0.984,1.00163,0.528269,0.517348,0.517987,0.431921,0.738932,0.271377];
tpRcdA =[0,9.4,26.7,44.8,62.1176,79.8118,97.5059,115.2,133.271,151.341]  ;

%PleD and PdeA
tpPdeA = [0 18.75 37.5 56.25 75 93.75 112.5 131.25 150];
dPdeA = [1  0.41  0.158  0.065  0.072  0.173  0.44   0.82  0.998]; %Abel 2011
%tpPdeA = [0 16.7 33.4 50.1 66.8 83.5 100.2 116.9 133.6 150]+1200;%Abel 2011
%dPdeA = [15600.731 3969.66 1534.024 1541.033 1240.74 1195.589 4955.116 6878.43 6766.468 19897.229];
tpPdeA = [15 30 45 60 75 90 105 120 135 150];
dPdeA = [0.741673635 0.316406924 0.107619862 0.046731771 0.046410539 0.114836436 0.348004248 0.626064795 0.798451444 1.00000012];

tpPleD=[0 10 20 30 40 50 60 70 80 90 100 110 120 130 140];
pPleD =[0.097 0.112 0.089 0.067 0.049 0.036 0.027 0.021 0.03 0.057 0.085 0.107 0.123 0.061 0.089];

tpcdG = [0    20    40    60    80   100   120 140]  ;
pcdG = [ 1/3   1   0.55  0.35   0.3  0.25  0.05  0.35]; %Abel 2013

tpcdG2 = [0    20    40    60    80   100   120]  ;
pcdG2 = [0.2292 1.0000 0.5000 0.2833 0.2250 0.1625 0];
% isolating the indices of the 7-10th cell cycle:
[~, a]=min(abs(tout(:)-1050));  % a and b are indeces of beginning and end of cell cycle 
[~, b]=min(abs(tout(:)-1500));

mut = 8
%% gcra mutant
if mut == 1  
    figure()
    f = gcf;    pbaspect([1.3 1 1])
    % set(f,'defaultAxesColorOrder',[left_color; right_color])
    box on;
    yyaxis left
    line( tout(a:b)-1050, yout(a:b, CtrA), 'Color', [0 .5 .5], 'LineWidth', 3, 'Linestyle', '-');
    line( tout(a:b)-1050, yout(a:b, DnaA), 'Color', [.5 .5 0], 'LineWidth', 3, 'Linestyle', '-');
    
    %hold on;
    ylabel('Normed. Conc. (CtrA, DnaA)', "FontSize",30)
    
    yyaxis right
    ylabel('Methylation Vars Count', "FontSize",30)
    line(tout(a:b) - 1050, yout(a:b, hcori), 'Color', 'k', 'LineWidth', 3, 'Linestyle', '-');
    line(tout(a:b) - 1050, yout(a:b, hCcrM), 'Color', 'm', 'LineWidth', 3, 'Linestyle', '-');
    line(tout(a:b) - 1050, yout(a:b, hCtrA), 'Color', 'b', 'LineWidth', 3, 'Linestyle', '--');
    xlim([0 450])
    ylim([0 2])
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    
    legend('CtrA','DnaA', 'h_{Cori}','h_{CcrM}','h_{ctrA}', 'best','FontSize',10)
    title('\it \DeltagcrA')
    xlabel('Time (min)','FontSize',24)
    
    set(f, 'MenuBar', 'figure');
end


%% ccrm mut
if mut == 2 
    figure()
    f = gcf;    pbaspect([1.3 1 1])
    % set(f,'defaultAxesColorOrder',[left_color; right_color])
    box on;
    yyaxis left
    line( tout(a:b)-1050, yout(a:b, DnaA), 'Color', [.5 .5 0], 'LineWidth', 3, 'Linestyle', '-');
    
    %hold on;
    ylabel('Normed. Conc. (DnaA)', "FontSize",30)
    
    yyaxis right
    ylabel('Methylation Vars Count', "FontSize",30)
    line(tout(a:b) - 1050, yout(a:b, hcori), 'Color', 'k', 'LineWidth', 3, 'Linestyle', '-');
    line(tout(a:b) - 1050, yout(a:b, hCcrM), 'Color', 'm', 'LineWidth', 3, 'Linestyle', '-');
    line(tout(a:b) - 1050, yout(a:b, hCtrA), 'Color', 'b', 'LineWidth', 3, 'Linestyle', '--');
    xlim([0 450])
    ylim([0 3])
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    
    legend('DnaA', 'h_{Cori}','h_{CcrM}','h_{ctrA}', 'best','FontSize',15)
    title('\it \DeltaccrM')
    xlabel('Time (min)','FontSize',24)
    
    set(f, 'MenuBar', 'figure');
end
%% cdG mut
if mut == 3 
    figure()
    f = gcf;    pbaspect([1.3 1 1])
    % set(f,'defaultAxesColorOrder',[left_color; right_color])
    box on;
    %yyaxis left
    line( tout(a:b)-1050, yout(a:b, CtrA), 'Color', 'k', 'LineWidth', 3, 'Linestyle', '-');
    line( tout(a:b)-1050, yout(a:b, CtrAP), 'Color', 'k', 'LineWidth', 3, 'Linestyle', '--');
    line( tout(a:b)-1050, yout(a:b, DnaA), 'Color', [.5 .5 0], 'LineWidth', 3, 'Linestyle', '-');
    
    hold on;
    ylabel('Normed. Conc.', "FontSize",30)
    
%     yyaxis right
%     ylabel('Methylation Vars Count', "FontSize",30)
%     line(tout(a:b) - 1050, yout(a:b, hcori), 'Color', 'k', 'LineWidth', 3, 'Linestyle', '-');
%     line(tout(a:b) - 1050, yout(a:b, hCcrM), 'Color', 'm', 'LineWidth', 3, 'Linestyle', '-');
%     line(tout(a:b) - 1050, yout(a:b, hCtrA), 'Color', 'b', 'LineWidth', 3, 'Linestyle', '--');
    xlim([0 450])
  %  ylim([0 3])
    ax = gca;
    ax.YAxis(1).Color = 'k';
   % ax.YAxis(2).Color = 'k';
    
    legend('CtrA', 'CtrA~P','DnaA', 'best','FontSize',15)
    title('\it cdG^0')
    xlabel('Time (min)','FontSize',24)
    
    set(f, 'MenuBar', 'figure');
end

%% pleD mut
if mut == 4 
    figure()
    f = gcf;    pbaspect([1.3 1 1])
    % set(f,'defaultAxesColorOrder',[left_color; right_color])
    box on;
    %yyaxis left
    line( tout(a:b)-1050, yout(a:b, CtrA), 'Color', [0 .5 .5], 'LineWidth', 3, 'Linestyle', '-');
    line( tout(a:b)-1050, yout(a:b, DnaA), 'Color', [.5 .5 0], 'LineWidth', 3, 'Linestyle', '-');
    
    hold on;
    ylabel('Normed. Conc.', "FontSize",30)
    
     %yyaxis right
     %ylabel("Normed. Conc. (cdG)")
     line(tout(a:b)-1050,yout(a:b,cdG), 'Color', [.5 0 .5], 'LineWidth', 3, 'Linestyle', '-');

    %ylim([0 .1])
    ax = gca;
    ax.YAxis(1).Color = 'k';
    %ax.YAxis(2).Color = 'k';
    
    legend('CtrA', 'DnaA','cdG', 'best','FontSize',15)
    title('\it \DeltapleD')
    xlabel('Time (min)','FontSize',24)
    
    set(f, 'MenuBar', 'figure');
end

%% pdeA mut
if mut == 5 
    figure()
    f = gcf;    pbaspect([1.3 1 1])
    % set(f,'defaultAxesColorOrder',[left_color; right_color])
    box on;
    %yyaxis left
    line( tout(a:b)-1050, yout(a:b, CtrA), 'Color', [0 .5 .5], 'LineWidth', 3, 'Linestyle', '-');
    
    hold on;
    ylabel('Normed. Conc.',"FontSize",30)
    
     %yyaxis right
     %ylabel("Normed. Conc. (DnaA, cdG)")
     line( tout(a:b)-1050, yout(a:b, DnaA), 'Color', [.5 .5 0], 'LineWidth', 3, 'Linestyle', '-');

     line(tout(a:b)-1050,yout(a:b,cdG), 'Color', [.5 0 .5], 'LineWidth', 3, 'Linestyle', '-');

    %ylim([0 4])
    ax = gca;
    ax.YAxis(1).Color = 'k';
   % ax.YAxis(2).Color = 'k';
    
    legend('CtrA', 'DnaA','cdG', 'best','FontSize',15)
    title('\it \DeltapdeA')
    xlabel('Time (min)','FontSize',24)
    
    set(f, 'MenuBar', 'figure');
end

%% dnaA mut
if mut == 6 
    figure()
    f = gcf;    pbaspect([1.3 1 1])
    % set(f,'defaultAxesColorOrder',[left_color; right_color])
    box on;
    yyaxis left
    line( tout(a:b)-1050, yout(a:b, DnaA), 'Color', [.5 .5 0], 'LineWidth', 3, 'Linestyle', '-');
    
    %hold on;
    ylabel('Normed. Conc. (DnaA)', "FontSize",30)
    
    yyaxis right
    ylabel('Methylation Vars Count', "FontSize",30)
    line(tout(a:b) - 1050, yout(a:b, hcori), 'Color', 'k', 'LineWidth', 3, 'Linestyle', '-');
    line(tout(a:b) - 1050, yout(a:b, hCcrM), 'Color', 'm', 'LineWidth', 3, 'Linestyle', '-');
    line(tout(a:b) - 1050, yout(a:b, hCtrA), 'Color', 'b', 'LineWidth', 3, 'Linestyle', '--');
    xlim([0 450])
    ylim([0 3])
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    
    legend('DnaA', 'h_{Cori}','h_{CcrM}','h_{ctrA}', 'best','FontSize',15)
    title('\it \DeltadnaA')
    xlabel('Time (min)','FontSize',24)
    
    set(f, 'MenuBar', 'figure');
end

%% SM921 mut
if mut == 7 
    figure()
    f = gcf;    pbaspect([1.3 1 1])
    % set(f,'defaultAxesColorOrder',[left_color; right_color])
    box on;
    yyaxis left
    line( tout(a:b)-1050, yout(a:b, CtrA), 'Color', [0 .5 .5], 'LineWidth', 3, 'Linestyle', '-');
    line( tout(a:b)-1050, yout(a:b, DnaA), 'Color', [.5 .5 0], 'LineWidth', 3, 'Linestyle', '-');
    
    %hold on;
    ylabel('Normed. Conc. (CtrA, DnaA)', "FontSize",30)
    
    yyaxis right
    ylabel('Methylation Vars Count', "FontSize",30)
    line(tout(a:b) - 1050, yout(a:b, hcori), 'Color', 'k', 'LineWidth', 3, 'Linestyle', '-');
    line(tout(a:b) - 1050, yout(a:b, hCcrM), 'Color', 'm', 'LineWidth', 3, 'Linestyle', '-');
    line(tout(a:b) - 1050, yout(a:b, hCtrA), 'Color', 'b', 'LineWidth', 3, 'Linestyle', '--');
    xlim([0 450])
    ylim([0 2])
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    
    legend('CtrA','DnaA', 'h_{Cori}','h_{CcrM}','h_{ctrA}', 'best','FontSize',10)
    title('SM921')
    xlabel('Time (min)','FontSize',24)
    
    set(f, 'MenuBar', 'figure');
end

%% ctra3 mut
if mut == 8 
    figure()
    f = gcf;    pbaspect([1.3 1 1])
    % set(f,'defaultAxesColorOrder',[left_color; right_color])
    box on;
    yyaxis left
    line( tout(a:b)-1050, yout(a:b, CtrA), 'Color', [0 .5 .5], 'LineWidth', 3, 'Linestyle', '-');
    line( tout(a:b)-1050, yout(a:b, DnaA), 'Color', [.5 .5 0], 'LineWidth', 3, 'Linestyle', '-');
    
    %hold on;
    ylabel('Normed. Conc. (CtrA, DnaA)', "FontSize",30)
    
    yyaxis right
    ylabel('Methylation Vars Count', "FontSize",30)
    line(tout(a:b) - 1050, yout(a:b, hcori), 'Color', 'k', 'LineWidth', 3, 'Linestyle', '-');
    line(tout(a:b) - 1050, yout(a:b, hCcrM), 'Color', 'm', 'LineWidth', 3, 'Linestyle', '-');
    line(tout(a:b) - 1050, yout(a:b, hCtrA), 'Color', 'b', 'LineWidth', 3, 'Linestyle', '--');
    xlim([0 450])
    ylim([0 2])
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    
    legend('CtrA','DnaA', 'h_{Cori}','h_{CcrM}','h_{ctrA}', 'best','FontSize',10)
    title('\it ctrA\Delta3\Omega')
    xlabel('Time (min)','FontSize',24)
    
    set(f, 'MenuBar', 'figure');
end
