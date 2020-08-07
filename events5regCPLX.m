%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% This is the script for detecting specific events of the cell cycle.
% This script is required by 'main5reg.m' for running.
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
function [value,isterminal,direction] = events5regCPLX(t, y)
% Locate the time when values passes through zero in a 
% increasing direction and stop integration.
global Ini Elong DNA Count hcori hCcrM hCtrA I mCcrM mDnaA ...
    mGcrA mSciP mCtrA CcrM DnaA GcrA SciP CtrA Sup;

value=[sign(y(Ini) -  0.05); sign(2*y(Elong) - 0.2*y(Count)); sign(2*y(Elong) - 0.375*y(Count)); ...
     sign(2*y(Elong) - y(Count)); sign(y(Sup) - 1.5); sign(y(Sup) - 0.3); sign(y(hcori)-0.05)];
isterminal=[1; 1; 1; 1; 1; 1; 1];
direction=[+1; +1; +1; +1; +1; +1; -1];
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
