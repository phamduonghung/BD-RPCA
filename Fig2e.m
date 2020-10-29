%%% code matlab of Fig2e-ground truth %%%%%                                              ;
clear  all;
close all 
%% Set Current Folder of MATLAB being BD-RPCA-GitHub and Add Path
addpath(genpath(fullfile(pwd)));
test=1;
nomfichier='simu_conv' 
result_folder = fullfile(pwd,'Results');
mkdir(result_folder)
%% Loading data
load_data_US;
%% ideal simulation figure
FigFeatures.title=1;
FigFeatures.result_folder = result_folder;
FigFeatures.mm=0;
FigFeatures.bar=1;
FigFeatures.print=0;
FigFeatures.nomtest = 'ground_truth';
IlogPWTD = -35*zeros(Nz,Nx,Nt);
IlogPWTD(200:270,45:58,:)=1;
IlogPWTD(145:180,90:100,:)=1;
FigFeatures.color1='b';
FigFeatures.color2= '--b';
%FigFeatures.bar=1;
[~,MIlogPWTD]=Dopplerplot(IlogPWTD,espace_xx,espace_zz,test,FigFeatures);                 
