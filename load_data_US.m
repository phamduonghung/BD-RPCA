%% Loading data US
iHS=0; % Not run Oleg
[~,name,~] = fileparts(nomfichier);
load(fullfile(pwd,'Data',sprintf('%s.mat',nomfichier)))                                ; %chargement de la matrice
if test ==1    
    [Nz,Nx,Nt] = size(M1)                           ; %Attribution de la taille de la matrice RF
    if 1
        for k = 1:Nt
            Moy = (mean(M1(:,:,k)))' * ones(1,Nz) ; 
            M1(:,:,k) = hilbert(M1(:,:,k)-Moy')              ; %application de la transform?e de Hilbert pour passer en donn?es complexes     
        end
        psf_fichier=fullfile(pwd,'Data','psf_simu.mat');
    end
    %%% INITIALISATION DE LA MATRICE PSF    
    if iHS==0
        psf_fichier=fullfile(pwd,'Data','psf_simu.mat');
    end
    load(psf_fichier);                                %chargement de la matrice                              
    % tic 
    [m_hh,n_hh] = size(psf); 
    shift_h = zeros(Nz*Nx, Nt);
    shift_h(1:m_hh,1:n_hh) = psf;
    H = fft2(circshift( shift_h, 1-[floor((m_hh+1)/2),floor((n_hh+1)/2)] ) );
    %     H = ones(Nz*Nx,Nt)+1i*ones(Nz*Nx,Nt);    
    espace_xx=1:Nx;
    espace_zz=1:Nz;
end