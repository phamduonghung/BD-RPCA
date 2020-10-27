%%% code matlab of Fig2b BD-RPCA%%%%%                                              ;
clear  all;
close all 
%% Let Current Folder of MATLAB contain BD-RPCA-GitHub and Add Path
addpath(genpath(fullfile(pwd)));

%% A modifier
test = 1; % For figure 2a of the paper, keep test=1
%%
if test ==1
    nomfichier='simu_conv' 
    seuil_tissu = 2;
    seuil_bruit = 15;
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
result_folder = fullfile(pwd,'Results');
mkdir(result_folder)
%% Loading data
iHS=0; % Not run Oleg
load_data_US;
[M,m,n,p] = convert_video3d_to_2d(M1);
%M = M/abs(max(M(:)));
%[DR,M]=dynarange(M); 
% show_2dvideo(M,m,n);

%%
tDRPCAStart = tic;           % pair 2: tic
fprintf('Running DRPCA....\n')
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

[T, S] = DRPCA(M,H,Lambda); %
Mfinale=reshape(S,Nz,Nx,Nt);
FigFeatures.nomtest = sprintf('DRPCA_%s',nomfichier); % Name 
Dopplerplot(Mfinale,espace_xx,espace_zz,test,FigFeatures); 
save(sprintf('%s/DRPCA_%s.mat', result_folder,nomfichier),'Mfinale')
clear Mfinale
tDRPCAEnd = toc(tDRPCAStart)      % pair 2: toc   