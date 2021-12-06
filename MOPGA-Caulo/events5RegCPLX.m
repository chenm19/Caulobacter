%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% This is the script for detecting specific events of the cell cycle.
% This script is required by 'main5RegCPLX.m' for running.
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% Locate the time when values passes through zero in a 
% increasing direction and stop integration.

function [value,isterminal,direction] = events5RegCPLX(t, y)
global Ini Elong Zring;

value=[sign(y(Ini) -  0.05); sign(y(Elong) - 0.2); sign(y(Elong) - 0.375); ...
     sign(y(Elong) - 1); sign(y(Zring) - 1); sign(y(Ini)-0.15)];
%  sign(y(Sup) - 0.3); 
%  sign(y(Ini)-0.15)
%  sign(y(hcori)-0.05)
isterminal=[1; 1; 1; 1; 1; 1;];
direction=[+1; +1; +1; +1; +1; +1;];
end
