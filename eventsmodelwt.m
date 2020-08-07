%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% Script for numericaly simulation of Caulobacter stalked cell cycle model  % 
% See the paper at the PLoS computational Biology' journal for details      % 
% Oct. 31, 2007. Shenghua Li, Paul Brazhnik, Bruno Sobral & John Tyson      %
%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
%
% This is the script for detecting specific events of stalked cell cycle.
% This script is required by 'caulomodelwt.m' for running.
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
function [value,isterminal,direction] = eventsmodewt(t, y)
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% Locate the time when values passes through zero in a increasing direction 
% and stop integration.
value=[sign(y(15)-0.05);sign(y(16)-0.2);sign(y(16)-0.375);sign(y(16)-0.625);sign(y(16)+1-y(18));sign(y(5)-0.9)];
isterminal=[1;1;1;1;1;1];
direction=[+1;+1;+1;+1;+1;+1];
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%

end
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%