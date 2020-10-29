%%% code matlab of Fig2b: BD-RPCA%%%%%                                              ;
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

%% Some figure parameters
FigFeatures.title=1;
FigFeatures.result_folder = result_folder;
FigFeatures.mm=0;
FigFeatures.bar=1;
FigFeatures.print=0;

%% APPLICATION DE LA BDRPCA
M = reshape(M1(:),Nz*Nx,Nt) ; %Construction de la matrice de Casorati

tBDRPCAStart = tic;           % pair 2: tic
%% Lambda1 Parameters
Lambda = 3./sqrt(max(Nz*Nx,Nt));
Lambda1 = 1./sqrt(max(Nz*Nx,Nt));
%%
tRPCAStart = tic;           % pair 2: tic
fprintf('Initialization RPCA....\n')
[T0, ~] = RobustPCA_Doppler(M,Lambda); %
tRPCAEnd = toc(tRPCAStart)      % pair 2: toc
%%
fprintf('Running estimated initial PSF ....\n')
max_iter = 3;
Mt = reshape(M-T0,Nz,Nx,Nt);
M11 = squeeze(mean(Mt,3));
[H,psf0] = Hestimate(M11,Nz,Nx,Nt);
fprintf('Initialized PSF size: %d-%d\n',size(psf0,1),size(psf0,2))
clear Mt M11 

%% Stop condition
tol  = 1e-3;
xtmp = M;
Ttmp = T0;
err = zeros(1,max_iter);
normM = norm(M, 'fro');

for iter = 1:max_iter
    iter
    fprintf('Running estimated DRPCA for iteration %d....\n',iter)
    [T, x] = DRPCA(M,H,Lambda1); % S <-> B (blood) and  L <->T (tissue) and M <-> S  and H<-> D in paper        
                   
    % Stop Condition
    Z1 = x-xtmp;    
    err(1,iter) = norm(Z1, 'fro') / normM; 
    xtmp=x;   
    
    if (err(1,iter) > tol)    
        Mt = reshape(M-T,Nz,Nx,Nt);
        M11 = squeeze(mean(Mt,3));
        fprintf('Running estimated PSF for iteration %d....\n',iter+1)
        [H,psf1] = Hestimate(M11,Nz,Nx,Nt); 
        fprintf('PSF size for iteration %d: %d-%d\n',iter+1,size(psf1,1),size(psf1,2))         
    else 
        break;
    end    
    clear Mt M11 psf1
end
tBDRPCAEnd = toc(tBDRPCAStart)      % pair 2: toc
%% AFFICHAGE DE L'IMAGE DEROULANTE SELON Nt APRES SEUILLAGE/FILTRAGE
Mfinale=reshape(x,Nz,Nx,Nt);
%save(sprintf('%s/BDRPCA_%s.mat', result_folder,nomfichier),'Mfinale')

%% Doppler de puissance
FigFeatures.nomtest = sprintf('BDRPCA_%s',nomfichier); % Name 
Dopplerplot(Mfinale,espace_xx,espace_zz,test,FigFeatures); 
clear Mfinale 

