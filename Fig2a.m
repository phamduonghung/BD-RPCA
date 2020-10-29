%%% code matlab of Fig2a: SVD%%%%%                                              ;
clear  all;
close all 
%% Set Current Folder of MATLAB being BD-RPCA-GitHub and Add Path
addpath(genpath(fullfile(pwd)));

%% A modifier
test = 1; % For figure 2a of the paper, keep test=1
%% Some parameters
nomfichier='simu_conv' 
seuil_tissu = 2;
seuil_bruit = 15;
result_folder = fullfile(pwd,'Results');
mkdir(result_folder)
%% Loading data
load_data_US;
[M,m,n,p] = convert_video3d_to_2d(M1);
%% SVD running
fprintf(sprintf('performing SVD...\n'))
tSVDStart = tic;           % pair 2: tic
Mnew = M'*M                 ; %Matrice carr?e
[V,D2,Vt] = svd(Mnew)       ; %Application de la SVD
D = sqrt(D2)                ; %Matrice des valeurs singuli?res
U = M*V/D                   ; %Calcul de la matrice spatiale des vecteurs singuliers
fprintf('Number of singular values: %d\n', length(diag(D)))

f=ones(1,Nt)                    ; %cr?ation d'un vecteur ones
f(1:seuil_tissu)=[0]            ; %Application du seuil tissu sur le vecteur 
f(seuil_bruit:Nt)=[0]           ; %Application du seuil bruit sur le vecteur
If=diag(f)                      ; %Matrice diagonale identit? filtr?e par les seuils
Mf=M*V*If*V'                    ; %Calcul de la matrice finale    

%% AFFICHAGE DE L'IMAGE DEROULANTE SELON Nt APRES SEUILLAGE/FILTRAGE
Mfinale=reshape(Mf,Nz,Nx,Nt)    ; 
%save(sprintf('%s/SVD_%s.mat', result_folder,nomfichier),'Mfinale')

%% Doppler de puissance
% Figures Parameters 
FigFeatures.title=1; % Figure title 0 ou 1
FigFeatures.result_folder = result_folder;
FigFeatures.mm=0; 
FigFeatures.bar=1; % Colorbar 0 or 1 
FigFeatures.print=0; % Pdf Figure Print: 0 or 1 through export_fig 
FigFeatures.nomtest = sprintf('SVD_%s',nomfichier); % Name 
Dopplerplot(Mfinale,espace_xx,espace_zz,test,FigFeatures); 
clear Mfinale 