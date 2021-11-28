%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% This is the script for reading file sOutput.txt generated from VTMOP.
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
close all; clear; clc;

format long
%% Input data.
disp('Begin to transfer ...')
fid = fopen('sOutput_T10k.txt');
%%%%% Skip 1 lines 
for i = 1:1
    tline = fgetl(fid);
    if ~ischar(tline)
        break;
    end
    disp(tline)
end

%%%%% Read the number of Row
tline = fgetl(fid);
RowNum = str2num(tline(46:end));

tline = fgetl(fid);
FunNum = str2num(tline(41:end));
%%%% Skip the following 3 lines
for i = 1:2
    tline = fgetl(fid);
    if ~ischar(tline)
        break;
    end
    disp(tline)
end

%%%% Read the Nondominated point set
Fval = fscanf(fid, '%f %f', [2, RowNum]);
Fval = Fval';

%%%% Skip the following 2 lines
for i = 1:3
    tline = fgetl(fid);
    if ~ischar(tline)
        break;
    end
    disp(tline)
end

% %%%%% Readout Efficient point set
% Val = fscanf(fid, '%f %f', [47, RowNum]);
% Val = Val';
% 
% %%%% Skip the following 3 lines
% for i = 1:1
%     tline = fgetl(fid);
%     if ~ischar(tline)
%         break;
%     end
%     disp(tline)
% end
% %%%%% Readout Full database of objective values observed
% Funval = fscanf(fid, '%f %f', [2, FunNum]);
% Funval = Funval';
% 
% %%%% Skip the following 3 lines
% for i = 1:3
%     tline = fgetl(fid);
%     if ~ischar(tline)
%         break;
%     end
%     disp(tline)
% end
% 
% %%%%% Readout Full database of evaluated design points
% FunctVal = fscanf(fid, '%f %f', [47, FunNum]);
% FunctVal = FunctVal';
% 
% fclose(fid);
% disp('It is the end of transfer ...')
figure(1)
plot(Fval(:,1),Fval(:,2),'o')

xlabel('Objective 1 - Exp(x)')
ylabel('Objective 2 - Exp(-x)')
% xlabel('Objective 1 - Protein level')
% ylabel('Objective 2 - Cell cycle')
title('Pareto front')
legend('VTMOP-VTDIRECT95')
