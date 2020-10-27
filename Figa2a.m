close all;
clear all;
clc;

%% Add Path
running_folder = 'C:\Users\dpham\ownCloud\Working\Atempo\';
addpath(genpath(fullfile(running_folder,'BD-RPCA-GitHub')));
addpath(genpath(fullfile(running_folder,'Data')));
%% A modifier
test = 3;
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
M = M/abs(max(M(:)));
%[DR,M]=dynarange(M); 
% show_2dvideo(M,m,n);

%%
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
% save(sprintf('%s/SVD.mat', result_folder),'Mfinale') % remove % if want to save the SVD result SVD.mat

%% Doppler de puissance
% Figures Parameters 
FigFeatures.title=1; % Figure title 0 ou 1
FigFeatures.result_folder = result_folder;
FigFeatures.mm=0; 
FigFeatures.bar=1; % Colorbar 0 or 1 
FigFeatures.print=0; % Figure 0 or 1
FigFeatures.nomtest = 'SVD'; % Name 
Dopplerplot(Mfinale,espace_xx,espace_zz,test,FigFeatures); 
clear Mfinale 