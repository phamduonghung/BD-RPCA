%%% code matlab of Fig2b: RPCA%%%%%                                              ;
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
iHS=0; % Not run Oleg
load_data_US;
[M,m,n,p] = convert_video3d_to_2d(M1);

%% RPCA running
tRPCAStart = tic;           % pair 2: tic
fprintf('Running RPCA....\n')
Lambda = 3./sqrt(max(Nz*Nx,Nt));
[T, S] = RobustPCA_Doppler(M,Lambda); %
tRPCAEnd = toc(tRPCAStart)      % pair 2: toc
%% AFFICHAGE DE L'IMAGE DEROULANTE SELON Nt APRES SEUILLAGE/FILTRAGE
Mfinale=reshape(S,Nz,Nx,Nt)    ; 
%save(sprintf('%s/RPCA_%s.mat', result_folder,nomfichier),'Mfinale')

%% Doppler de puissance
% Figures Parameters 
FigFeatures.title=1; % Figure title 0 ou 1
FigFeatures.result_folder = result_folder;
FigFeatures.mm=0; 
FigFeatures.bar=1; % Colorbar 0 or 1 
FigFeatures.print=0; % Pdf Figure Print: 0 or 1 through export_fig 
FigFeatures.nomtest = sprintf('RPCA_%s',nomfichier); % Name 
Dopplerplot(Mfinale,espace_xx,espace_zz,test,FigFeatures); 
clear Mfinale 