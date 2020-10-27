close all;
clear all;
clc;

%% Add Path
running_folder = 'C:\Users\dpham\ownCloud\Working\Atempo\';
addpath(genpath(fullfile(running_folder,'BD-RPCA-GitHub')));
addpath(genpath(fullfile(running_folder,'Data')));
%% A modifier
test = 1; % For figure 2a of the paper, keep test==1
mm=0;
%%
if test ==1
    nomfichier='simu_conv' 
    seuil_tissu = 2;
    seuil_bruit = 15;
    mm=1;
elseif test ==2
    nomfichier='cerveau_sain'
    seuil_tissu = 100;
    seuil_bruit = 150;
elseif test ==3
    nomfichier='peri' 
    seuil_tissu = 100;
    seuil_bruit = 150;
else 
    nomfichier='tumeur'
    seuil_tissu = 100;
    seuil_bruit = 200;
end
result_folder = fullfile(running_folder,'BD-RPCA-GitHub','Results');
mkdir(result_folder)
%% Loading data
iHS=0; % Not run Oleg
load_data_US;
[M,m,n,p] = convert_video3d_to_2d(M1);
%M = M/abs(max(M(:)));
%[DR,M]=dynarange(M); 
% show_2dvideo(M,m,n);

%%
tRPCAStart = tic;           % pair 2: tic
fprintf('Running RPCA....\n')
if test==1
    Lambda = 3./sqrt(max(Nz*Nx,Nt));
    Lambda1 = 1./sqrt(max(Nz*Nx,Nt));
elseif test==2
    Lambda = 1.3*1./sqrt(max(Nz*Nx,Nt));
    Lambda1 = 1*1./sqrt(max(Nz*Nx,Nt)); % with C=0.5*[1 1]; B=[32 24];                                  
elseif test==3
    Lambda = 0.9*1./sqrt(max(Nz*Nx,Nt));
    Lambda1 = 1.15./sqrt(max(Nz*Nx,Nt));
else         
    Lambda = 1.2*1./sqrt(max(Nz*Nx,Nt));
    Lambda1 = 1.1*1./sqrt(max(Nz*Nx,Nt)); % with C=0.5*[1 1]; B=[32 24];  
end

[TtrueR, StrueR] = RobustPCA_Doppler(M,Lambda); %
tRPCAEnd = toc(tRPCAStart)      % pair 2: toc
%% AFFICHAGE DE L'IMAGE DEROULANTE SELON Nt APRES SEUILLAGE/FILTRAGE
Mfinale=reshape(StrueR,Nz,Nx,Nt)    ; 
%save(sprintf('%s/RPCA_simu.mat', result_folder),'Mfinale')

%% Doppler de puissance
% Figures Parameters 
FigFeatures.title=1; % Figure title 0 ou 1
FigFeatures.result_folder = result_folder;
FigFeatures.mm=0; 
FigFeatures.bar=1; % Colorbar 0 or 1 
FigFeatures.print=1; % Pdf Figure Print 0 or 1
FigFeatures.nomtest = 'RPCA_simu'; % Name 
Dopplerplot(Mfinale,espace_xx,espace_zz,test,FigFeatures); 
clear Mfinale 