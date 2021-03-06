%%% code matlab of Fig2b: DRPCA%%%%%                                              ;
clear  all;
close all 
%% Set Current Folder of MATLAB being BD-RPCA-GitHub and Add Path
addpath(genpath(fullfile(pwd)));

%% Some parameters
test = 1; % For figure 2a of the paper, keep test=1 
nomfichier='simu_conv' 
result_folder = fullfile(pwd,'Results');
mkdir(result_folder)
%% Loading data
load_data_US;
[M,m,n,p] = convert_video3d_to_2d(M1);

%% BDRPCA running
tDRPCAStart = tic;           % pair 2: tic
fprintf('Running DRPCA....\n')
Lambda = 3./sqrt(max(Nz*Nx,Nt));
[T, S] = DRPCA(M,H,Lambda); %
Mfinale=reshape(S,Nz,Nx,Nt);
tDRPCAEnd = toc(tDRPCAStart)     
%save(sprintf('%s/DRPCA_%s.mat', result_folder,nomfichier),'Mfinale')

%% Figures Parameters 
FigFeatures.title=1; % Figure title 0 ou 1
FigFeatures.result_folder = result_folder;
FigFeatures.mm=0; 
FigFeatures.bar=1; % Colorbar 0 or 1 
FigFeatures.print=0; % Pdf Figure Print 0 or 1 through export_fig 
FigFeatures.nomtest = sprintf('DRPCA_%s',nomfichier); % Name 
Dopplerplot(Mfinale,espace_xx,espace_zz,test,FigFeatures); 
clear Mfinale
