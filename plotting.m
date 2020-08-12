%% Ploting the time courses of model variables
close all; clf;

load('output.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Five transcriptional regulators
% figure(11)
% title('Four transcriptional regulators from Shenghua model')
% p1 = line(tout, yout(:, CtrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
% p2 = line(tout, yout(:, DnaA), 'Color', 'm', 'LineWidth', 2, 'Linestyle', '-');
% p3 = line(tout, yout(:, GcrA), 'Color', 'b', 'LineWidth', 2, 'Linestyle', '-');
% p4 = line(tout, yout(:, CcrM), 'Color', 'y', 'LineWidth', 2, 'Linestyle', '-');
% p5 = line(tout, yout(:, SciP), 'Color', 'r', 'LineWidth', 2, 'Linestyle', '-');
% h = legend('CtrA', 'DnaA',  'GcrA', 'CcrM', 'SciP', 'Location', 'North');

% figure(22)
% title('Five transcriptional regulators from Shenghua model')
% p1 = line(tout, yout(:, mCtrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
% p2 = line(tout, yout(:, mDnaA), 'Color', 'm', 'LineWidth', 2, 'Linestyle', '-');
% p3 = line(tout, yout(:, mGcrA), 'Color', 'b', 'LineWidth', 2, 'Linestyle', '-');
% p4 = line(tout, yout(:, mCcrM), 'Color', 'y', 'LineWidth', 2, 'Linestyle', '-');
% p5 = line(tout, yout(:, mSciP), 'Color', 'g', 'LineWidth', 2, 'Linestyle', '-');
% h = legend('mCtrA', 'mDnaA',  'mGcrA', 'mCcrM', 'mSciP','Location', 'North');

% hold on;

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
tot = 28;

CPLX1 = 23;
CpdR = 24;
CpdRP = 25;
CPLX2 = 26;
RcdA = 27;
CPLX3 = 28;

t =  [5, 30, 60 90 120 150] + 1050;
dDnaA = [1262	1164 500	496 	1066	1028]./445;
dCcrM = [17	17	24	445	435	101]./445;
dSciP = [2858	451	199	796	2956	2649]./445;
dGcrA = [550	3001	1275	476	725	1493]./445;
dCtrA = [195	612	2251	3216	2486	1617]./445;


tp = [0, 20, 40, 60, 80, 100, 120, 140] + 1050;
tpCcrM = [0, 20, 40, 60, 80, 100, 120] + 1050;
tp2 = [9, 27, 44, 62, 80, 98, 117, 134, 152] + 1050;

pDnaA = 2*[0.4567	0.9857	0.9147	0.3154	0.128	0.0793	0.4203	0.6582];
pGcrA = 6*[0	0.5021	1	0.884	1	0.554	0.662	0.7388];
pCtrA1 = 7*[0.7791	0	0	0	0.4682	1	0.8591	0.6672];
pCcrM = [0.979012341	0.245634747	0.105451343	0.054014407	0.031971569	0.078998499	1];
pSciP = 0.7*[7.757911401	7.931725156	4.139518022	1	0.446070124	0.429447455	0.761125117	3.856508369	8.881823903];
pCtrA2 = [0.8158	0.614	0.1881	0.0323	0.3822	0.7033	0.8924	1	0.782];

%% figure 1: plotting concentration of 5 main regulators mRNA against time
figure(1)

% ccrM mRNA subplot
subplot(2,3,1)
p1 = line(tout, yout(:, mCcrM), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
plot(t, dCcrM, 'r*-')
legend('Simulation','Experiment')
title('\it{ccrM}')

% dnaA mRNA subplot
subplot(2,3,2)
p1 = line(tout, yout(:, mDnaA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
plot(t, dDnaA, 'r*-')
legend('Simulation','Experiment')
title('\it{dnaA}')

% gcra mRNA subplot
subplot(2,3,3)
p1 = line(tout, yout(:, mGcrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
plot(t, dGcrA, 'r*-')
legend('Simulation','Experiment')
title('\it{gcrA}')

% sciP mRNA subplot
subplot(2,3,4)
p1 = line(tout, yout(:, mSciP), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
plot(t, dSciP, 'r*-')
legend('Simulation','Experiment')
title('\it{sciP}')

% ctra mRNA subplot
subplot(2,3,5)
p1 = line(tout, yout(:, mCtrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
plot(t, dCtrA, 'r*-')
legend('Simulation','Experiment')
title('\it{ctrA}')

%% figure 2: plotting concentration of 5 main regulator proteins (and protease complex) against time
figure(2)

% Ccrm protein subplot
subplot(2,3,1)
p1 = line(tout, yout(:, CcrM), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
plot(tpCcrM, pCcrM, 'r*-')
legend('Simulation','Experiment')
title('CcrM')

% DnaA protein subplot
subplot(2,3,2)
p1 = line(tout, yout(:, DnaA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
plot(tp, pDnaA, 'r*-')
legend('Simulation','Experiment')
title('DnaA')

% GcrA protein subplot
subplot(2,3,3)
p1 = line(tout, yout(:, GcrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
plot(tp, pGcrA, 'r*-')
legend('Simulation','Experiment')
title('GcrA')

% SciP protein subplot
subplot(2,3,4)
p1 = line(tout, yout(:, SciP), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
plot(tp2, pSciP, 'r*-')
legend('Simulation','Experiment')
title('SciP')

% CtrA protein subplot
subplot(2,3,5)
p1 = line(tout, yout(:, CtrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
plot(tp, pCtrA1, 'r*-')
legend('Simulation','Experiment')
title('CtrA')

% Complex3 protease subplot
subplot(2,3,6)
p1 = line(tout, yout(:, RcdA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
% plot(tp, pCtrA1, 'r*-')
legend('Simulation')
title('CPLX3')

% scatter(t, dCtrA, 'k*');
% scatter(t, dDnaA, 'm*');
% scatter(t, dGcrA, 'b*');
% scatter(t, dCcrM, 'y*');
% scatter(t, dSciP, 'g*');
% legend('Simulation','Experiment')
% title('mRNA of CtrA')
% axis([0 450 0 8])

%% Plotting DNA synthesis and Methylation variables

figure(3)

subplot(1,2,1)
% title('Four transcriptional regulators from Shenghua model')
p1 = line(tout, yout(:, hcori), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
p2 = line(tout, yout(:, hCcrM), 'Color', 'm', 'LineWidth', 2, 'Linestyle', '-');
p3 = line(tout, yout(:, hCtrA), 'Color', 'b', 'LineWidth', 2, 'Linestyle', '--');
h = legend('h_{Cori}', 'h_{ccrM}',  'h_{ctrA}', 'Location', 'North');
xlabel('time (min)')
% ylabel('Count')
% axis([0 450 0 1.5])

subplot(1,2,2)
% figure(44)
p2 = line(tout, yout(:, Elong), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
p3 = line(tout, yout(:, DNA), 'Color', 'm', 'LineWidth', 2, 'Linestyle', '-');
p4 = line(tout, yout(:, Count), 'Color', 'b', 'LineWidth', 2, 'Linestyle', '--');
% p5 = line(tout, yout(:, Ini), 'Color', 'r', 'LineWidth', 2, 'Linestyle', '-');
h = legend( 'Elongation',  'DNA', 'Chromosome');
% axis([0 600 0 0.1])

% figure(55)
% p1 = line(tout, yout(:, Elong), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
% axis([0 420 0 0.1])

% figure(66)
% % subplot(2,3,1)
% figure(61)
% avg = mean(yout(:, CcrM))*ones(1,30);
% space1 = zeros(1,90);
% space2 = zeros(1,30);
% et = 1:1:150;
% e1 = [space1,  avg,  space2];
% scatter(et, e1, 'r*', 'LineWidth', 4);
% hold on;
% ybar = zeros(size(tout));
% ybar(yout(:, CcrM)>mean(yout(:, CcrM)))=mean(yout(:, CcrM));
% scatter(tout, ybar, '.b', 'LineWidth', 4);
% plot(tout, yout(:,CcrM), 'k')
% title('CcrM')
% legend('Experiment','Simulation','Curve')
% axis([0 150 0 1])
% grid on 
% set(gca,'xtick',[0:30:150])
% 
% 
% % subplot(2,3,2)
% figure(62)
% avg = mean(yout(:, DnaA))*ones(1,50);
% space1 = zeros(1,10);
% space2 = zeros(1,90);
% e1 = [space1,  avg,  space2];
% scatter(et, e1, 'r*', 'LineWidth', 4);
% hold on;
% ybar = zeros(size(tout));
% ybar(yout(:, DnaA)>mean(yout(:, DnaA)))=mean(yout(:, DnaA));
% p1 = scatter(tout, ybar, 'b.', 'LineWidth', 4);
% plot(tout, yout(:,DnaA), 'k')
% title('DnaA')
% legend('Experiment','Simulation','Curve')
% axis([0 150 0 3])
% grid on 
% set(gca,'xtick',[0:30:150])
% 
% % subplot(2,3,3)
% figure(63)
% avg = mean(yout(:, GcrA))*ones(1,55);
% space1 = zeros(1,45);
% space2 = zeros(1,50);
% e1 = [space1,  avg,  space2];
% scatter(et, e1, 'r*', 'LineWidth', 4);
% hold on;
% ybar = zeros(size(tout));
% ybar(yout(:, GcrA)>mean(yout(:, GcrA)))=mean(yout(:, GcrA));
% p1 = scatter(tout, ybar, 'b.', 'LineWidth', 4);
% plot(tout, yout(:,GcrA), 'k')
% title('GcrA')
% legend('Experiment','Simulation','Curve')
% axis([0 150 0 6])
% grid on 
% set(gca,'xtick',[0:30:150])
% 
% % subplot(2,3,4)
% figure(64)
% avg = mean(yout(:, SciP))*ones(1,15);
% avg2 = mean(yout(:, SciP))*ones(1,60);
% space1 = zeros(1,85);
% e1 = [avg, space1, avg2];
% scatter(et, e1, 'r*', 'LineWidth', 4);
% hold on;
% ybar = zeros(size(tout));
% ybar(yout(:, SciP)>mean(yout(:, SciP)))=mean(yout(:, SciP));
% p1 = scatter(tout, ybar, 'b.', 'LineWidth', 4);
% plot(tout, yout(:,SciP), 'k')
% title('SciP')
% legend('Experiment','Simulation','Curve')
% axis([0 150 0 7])
% grid on 
% set(gca,'xtick',[0:30:150])
% 
% % subplot(2,3,5)
% figure(65)
% avg = mean(yout(:, CtrA))*ones(1,15);
% avg2 = mean(yout(:, CtrA))*ones(1,60);
% space1 = zeros(1,75);
% e1 = [avg, space1, avg2];
% scatter(et, e1, 'r*', 'LineWidth', 4);
% hold on;
% ybar = zeros(size(tout));
% ybar(yout(:, CtrA)>mean(yout(:, CtrA)))=mean(yout(:, CtrA));
% p1 = scatter(tout, ybar, 'b.', 'LineWidth', 4);
% plot(tout, yout(:,CtrA), 'k')
% title('CtrA')
% legend('Experiment','Simulation','Curve')
% axis([0 150 0 8])

% grid on 
% set(gca,'xtick',[0:30:150])


% h = legend('CtrA', 'DnaA',  'GcrA', 'CcrM', 'Location', 'North');
% axis([0 500 0.2 4])

% set(gca,'YTick',[1:4]) 
% lab = ['CcrM';'GcrA'; 'DnaA';'CtrA'];
% set(gca,'yticklabel',lab);
% title('Temporal control results from Shenghua model')
%%


%end of plotting
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%

% figure(66)
% subplot(5,1,1)
% % figure(61)
% avg = ones(1,30);
% space1 = zeros(1,90);
% space2 = zeros(1,30);
% et = 1:1:150;
% e1 = [space1,  avg,  space2];
% scatter(et, e1, 'or', 'LineWidth', 1, 'MarkerFaceColor', 'r' );
% hold on;
% ybar = zeros(size(tout));
% ybar(yout(:, CcrM)>mean(yout(:, CcrM)))=2;
% scatter(tout, ybar, 'sb', 'LineWidth', 2, 'MarkerFaceColor', 'b' );
% title('CcrM')
% % legend('Experiment','Simulation','Curve')
% axis([0 150 0.5 2.5])
% grid on 
% set(gca,'xtick',[0:30:150])
% set(gca,'ytick',[1,2])
% set(gca,'yticklabel',{'Experiment','Simulation'})
% 
% 
% subplot(5,1,2)
% % figure(62)
% avg = ones(1,50);
% space1 = zeros(1,10);
% space2 = zeros(1,90);
% e1 = [space1,  avg,  space2];
% scatter(et, e1, 'or', 'LineWidth', 1, 'MarkerFaceColor', 'r' );
% hold on;
% ybar = zeros(size(tout));
% ybar(yout(:, DnaA)>mean(yout(:, DnaA)))=2;
% scatter(tout, ybar, 'sb', 'LineWidth', 2, 'MarkerFaceColor', 'b' );
% title('DnaA')
% axis([0 150 0.5 2.5])
% grid on 
% set(gca,'xtick',[0:30:150])
% set(gca,'ytick',[1,2])
% set(gca,'yticklabel',{'Experiment','Simulation'})
% 
% subplot(5,1,3)
% % figure(63)
% avg = ones(1,55);
% space1 = zeros(1,45);
% space2 = zeros(1,50);
% e1 = [space1,  avg,  space2];
% scatter(et, e1, 'or', 'LineWidth', 1, 'MarkerFaceColor', 'r' );
% hold on;
% ybar = zeros(size(tout));
% ybar(yout(:, GcrA)>mean(yout(:, GcrA)))=2;
% scatter(tout, ybar, 'sb', 'LineWidth', 2, 'MarkerFaceColor', 'b' );
% title('GcrA')
% % legend('Experiment','Simulation','Curve')
% axis([0 150 0.5 2.5])
% grid on 
% set(gca,'xtick',[0:30:150])
% set(gca,'ytick',[1,2])
% set(gca,'yticklabel',{'Experiment','Simulation'})
% 
% subplot(5,1,4)
% % figure(64)
% avg = ones(1,15);
% avg2 = ones(1,60);
% space1 = zeros(1,75);
% e1 = [avg, space1, avg2];
% scatter(et, e1, 'or', 'LineWidth', 1, 'MarkerFaceColor', 'r' );
% hold on;
% ybar = zeros(size(tout));
% ybar(yout(:, SciP)>mean(yout(:, SciP)))=2;
% scatter(tout, ybar, 'sb', 'LineWidth', 2, 'MarkerFaceColor', 'b' );
% title('SciP')
% % legend('Experiment','Simulation','Curve')
% axis([0 150 0.5 2.5])
% grid on 
% set(gca,'xtick',[0:30:150])
% set(gca,'ytick',[1,2])
% set(gca,'yticklabel',{'Experiment','Simulation'})
% 
% subplot(5,1,5)
% % figure(65)
% avg = ones(1,15);
% avg2 = ones(1,90);
% space1 = zeros(1,45);
% e1 = [avg, space1, avg2];
% scatter(et, e1, 'or', 'LineWidth', 1, 'MarkerFaceColor', 'r' );
% hold on;
% ybar = zeros(size(tout));
% ybar(yout(:, CtrA)>mean(yout(:, CtrA)))=2;
% scatter(tout, ybar, 'sb', 'LineWidth', 2, 'MarkerFaceColor', 'b' );
% title('CtrA')
% % legend('Experiment','Simulation','Curve')
% axis([0 150 0.5 2.5])
% grid on 
% set(gca,'xtick',[0:30:150])
% set(gca,'ytick',[1,2])
% set(gca,'yticklabel',{'Experiment','Simulation'})


%%  Complex figures 
% plot Complex1, CpdR, CpdRP vs Time
figure();
hold on; 
subplot(3,1,1);
plot(tout, yout(:, CPLX1),'k');
xlabel('Time/min')
ylabel('Complex1')
subplot(3,1,2);
plot(tout,yout(:,CpdR),'r');
xlabel('Time/min')
ylabel('CpdR')
subplot(3,1,3);
plot(tout,yout(:,CpdRP),'b');
xlabel('Time/min')
ylabel('CpdRP')

% plot Complex2, RcdA vs Time
figure();
hold on;
subplot(2,1,1);
plot(tout,yout(:,CPLX2),'k');
xlabel('Time/min')
ylabel('Complex2')
subplot(2,1,2);
plot(tout,yout(:,RcdA),'r');
xlabel('Time/min')
ylabel('RcdA')

% plot cpdR + cpdRP vs Time
figure();
hold on;
plot(tout,yout(:,CpdR)+yout(:,CpdRP))
xlabel('Time/min')
ylabel('total CpdR')

% plot Complex3 vs Time
figure();
hold on;
plot(tout,yout(:,CPLX3))
xlabel('Time/min')
ylabel('Complex3')

%% Comparing experimental vs simulated Complex data
% plotting the 4th cell cycle, 450-600 minutes

figure();
[~, a]=min(abs(tout(:)-450));  % a and b are indeces of beginning and end of cell cycle 
[~, b]=min(abs(tout(:)-600));

% subplot(3,1,1);
%time=[0,20,40,60,80,100,120];
% % cpdrp=[9118.983,4425.406,3905.163,1912.92,5768.719,8952.497,6521.347];
% cpdrp=[9118.983;4425.406;3905.163;1912.92;5768.719;8952.497;22145.296];
% % cpdrp=cpdrp/max(cpdrp);
% dif=max(cpdrp)-min(cpdrp);
% cpdrp=(cpdrp-min(cpdrp))/dif;
% 
% scatter(time+20,cpdrp,'ro','MarkerFaceColor','r')
% hold on;
% CPDRP=yout(:,3);
% % plot(tout,CPDRP);
% CPDRP=CPDRP(a:b);
% DIF=max(CPDRP)-min(CPDRP);
% CPDRP=(CPDRP-min(CPDRP))/DIF;
% plot(tout(a:b)-450,CPDRP)%plot the forth cell cycle 450-600min
% xlabel('Time/min')
% ylabel('CpdRP')
% legend('experimental data','simulated CpdRP')
% hold on;


subplot(2,1,1)
% cpdr=[21001.196,20613.125,13581.933,10400.397,10563.811,13216.569,20216.276];
cpdr=[27532.359,40939.622,20027.844,10400.397,10563.811,13216.569,20216.276];
% cpdr=cpdr/max(cpdr);
dif=max(cpdr) - min(cpdr);
cpdr=(cpdr-min(cpdr))/dif;
scatter(time+10,cpdr,'ro','MarkerFaceColor','r')  % plotting experimental cpdr points
hold on;
CPDR = yout(:, CpdR);                                   
CPDR = CPDR(a:b);                  % gathering relevant simulated CpdR data
% CPDR=CPDR/max(CPDR);
DIF = max(CPDR) - min(CPDR);
CPDR = (CPDR-min(CPDR))/DIF;
plot(tout(a:b)-450, CPDR)            % plotting simulated CpdR
xlabel('Time/min')
ylabel('CpdR')
legend('experimental data', 'simulated CpdR', 'location', 'northeastoutside')
hold on;

time2=[0,9.4,26.7,44.8,62.1176,79.8118,97.5059,115.2,133.271,151.341];
rcda=[0.0000136,0.093,0.984,1.00163,0.528269,0.517348,0.517987,0.431921,0.738932,0.271377];
% rcda=rcda/max(rcda);
dif=max(rcda)-min(rcda);
rcda=(rcda-min(rcda))/dif;

subplot(2, 1, 2)
scatter(time2, rcda, 'ro', 'MarkerFaceColor', 'r')  % plotting experimental rcda points
hold on;
RCDA = yout(:, RcdA);
RCDA = RCDA(a:b);                    % gathering relevant simulated CpdR data
% RCDA=RCDA/max(RCDA);
DIF = max(RCDA) - min(RCDA);
RCDA = (RCDA-min(RCDA))/DIF;         % plotting simulated RcdA data
plot(tout(a:b)-450, RCDA)
xlabel('Time/min')
ylabel('RcdA')
legend('experimental data','simulated RcdA', 'location', 'northeastoutside')


