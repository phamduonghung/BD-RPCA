%% Loading data US
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
else 
    M1=real(double(IQ));
    espace_xx=espace_x;
    espace_zz=espace_z;
    % Attribution de la taille de la matrice RF
    [Nz,Nx,Nt] = size(M1); 
    %%% INITIALISATION DE LA MATRICE PSF
    if test ==2 
        test
        psf_fichier=fullfile(pwd,'Data','psf_cerveau.mat');
        load(psf_fichier);                                %chargement de la matrice                             
        psf = real(double(psf(:,:)));       
        [m_hh,n_hh] = size(psf); 
        shift_h = zeros(Nz*Nx, Nt);
        shift_h(1:m_hh,1:n_hh) = psf;
        H = fft2(circshift(shift_h, 1-[floor((m_hh+1)/2),floor((n_hh+1)/2)]));          
    elseif test ==3
        psf_fichier=fullfile(pwd,'Data','psf_peri.mat');
        load(psf_fichier);                                %chargement de la matrice                             
        psf = real(double(psf(:,:)));
        [m_hh,n_hh] = size(psf); 
        shift_h = zeros(Nz*Nx, Nt);
        shift_h(1:m_hh,1:n_hh) = psf;
        H = fft2( circshift( shift_h, 1-[floor((m_hh+1)/2),floor((n_hh+1)/2)] ) );
    else
        psf_fichier=fullfile(pwd,'Data','psf_tumeur.mat');
        load(psf_fichier);                                %chargement de la matrice                             
        apsf = real(double(IQ(:,:,1)));
        [Pz,Px,Pt] = size(apsf); 
        % choisir le point du milieu
        m_h=32;
        n_h=16;
        point=[floor(Pz/2)-42 , floor(Px/2)+05];
        psf = apsf(point(1)-floor(m_h/2):point(1)+floor(m_h/2) , point(2)-floor(n_h/2):point(2)+floor(n_h/2));
        psf = psf/sum(psf(:));  
        %    tic 
        [m_hh,n_hh] = size(psf); 
        shift_h = zeros(Nz*Nx, Nt);
        shift_h(1:m_hh,1:n_hh) = psf;
        H = fft2( circshift( shift_h, 1-[floor((m_hh+1)/2),floor((n_hh+1)/2)] ) );
    end
end
