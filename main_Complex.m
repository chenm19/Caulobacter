%temperal_odes.m
clc;
clear all;
close all;

global p;

param_Complex();

% [T,Y]=ode15s('odes_Complex',[0,300],[0,1,1]);%Complex1, CpdR, CpdRP;
% [T,Y]=ode15s('odes_Complex',[0,900],[0.1,4,0,0,2]);%Complex2, RcdA
[T,Y]=ode15s('odes_Complex',[0,1500],[0.1,3,1,0,2,0]);%Complex3

figure();
hold on; 
subplot(3,1,1);
plot(T,Y(:,1),'k');
xlabel('Time/min')
ylabel('Complex1')
subplot(3,1,2);
plot(T,Y(:,2),'r');
xlabel('Time/min')
ylabel('CpdR')
subplot(3,1,3);
plot(T,Y(:,3),'b');
xlabel('Time/min')
ylabel('CpdRP')

figure();
hold on;
subplot(2,1,1);
plot(T,Y(:,4),'k');
xlabel('Time/min')
ylabel('Complex2')
subplot(2,1,2);
plot(T,Y(:,5),'r');
xlabel('Time/min')
ylabel('RcdA')

figure();
hold on;
plot(T,Y(:,2)+Y(:,3))
xlabel('Time/min')
ylabel('total CpdR')

figure();
hold on;
plot(T,Y(:,6))
xlabel('Time/min')
ylabel('Complex3')

%% COMPARISON
figure();
[~, a]=min(abs(T(:)-450));
[~, b]=min(abs(T(:)-600));
hold on;
% subplot(3,1,1);
time=[0,20,40,60,80,100,120];
% % cpdrp=[9118.983,4425.406,3905.163,1912.92,5768.719,8952.497,6521.347];
% cpdrp=[9118.983;4425.406;3905.163;1912.92;5768.719;8952.497;22145.296];
% % cpdrp=cpdrp/max(cpdrp);
% dif=max(cpdrp)-min(cpdrp);
% cpdrp=(cpdrp-min(cpdrp))/dif;
% 
% scatter(time+20,cpdrp,'ro','MarkerFaceColor','r')
% hold on;
% CPDRP=Y(:,3);
% % plot(T,CPDRP);
% CPDRP=CPDRP(a:b);
% DIF=max(CPDRP)-min(CPDRP);
% CPDRP=(CPDRP-min(CPDRP))/DIF;
% plot(T(a:b)-450,CPDRP)%plot the forth cell cycle 450-600min
% xlabel('Time/min')
% ylabel('CpdRP')
% legend('experimental data','simulated CpdRP')
% hold on;


subplot(2,1,1)
% cpdr=[21001.196,20613.125,13581.933,10400.397,10563.811,13216.569,20216.276];
cpdr=[27532.359,40939.622,20027.844,10400.397,10563.811,13216.569,20216.276];
% cpdr=cpdr/max(cpdr);
dif=max(cpdr)-min(cpdr);
cpdr=(cpdr-min(cpdr))/dif;

scatter(time+10,cpdr,'ro','MarkerFaceColor','r')
hold on;
CPDR=Y(:,2);
CPDR=CPDR(a:b);
% CPDR=CPDR/max(CPDR);
DIF=max(CPDR)-min(CPDR);
CPDR=(CPDR-min(CPDR))/DIF;

plot(T(a:b)-450,CPDR)
xlabel('Time/min')
ylabel('CpdR')
legend('experimental data','simulated CpdR')
hold on;
time2=[0,9.4,26.7,44.8,62.1176,79.8118,97.5059,115.2,133.271,151.341];
rcda=[0.0000136,0.093,0.984,1.00163,0.528269,0.517348,0.517987,0.431921,0.738932,0.271377];
% rcda=rcda/max(rcda);
dif=max(rcda)-min(rcda);
rcda=(rcda-min(rcda))/dif;

subplot(2,1,2)
scatter(time2,rcda,'ro','MarkerFaceColor','r')
hold on;
RCDA=Y(:,5);
RCDA=RCDA(a:b);
% RCDA=RCDA/max(RCDA);
DIF=max(RCDA)-min(RCDA);
RCDA=(RCDA-min(RCDA))/DIF;
plot(T(a:b)-450,RCDA)

xlabel('Time/min')
ylabel('RcdA')
legend('experimental data','simulated RcdA')
