function [H,psf] = Hestimate(X,Nz,Nx,Nt)
    %clc
    %clear
    %close all

    %------------------------- ESTIMATION PARAMETERS -------------------------
    %X=real(M11);
    C=0.99*[1 1];                               % LPF cut-off frequencies
    B=[64 48];                                    % block size for PSE
    L=[4 4];                                     % #splines for PSE

    lmd=2e-3;                                   % regularization parameter
    a=0.05;                                     % Huber parameter

    %-------------------- DATA READING AND NORMALIZATION  --------------------
    %d.rf = X;
    %d=load('data.mat');                         % reading RF data
    [DR,X]=dynarange(X);                   % normalization of RF data

    %----------------------------- DEMODULATION ------------------------------
    f0 = 5.625*10^6;
    fs = 25*10^6;
    fn=f0./fs;                              % normalized central frequency
    [iq,D]=rf2iq(X,fn,C);                    % demodulation of RF data
    
%     [iq1,D1]=rf2iq(real(X),fn,C);                    % demodulation of RF data
%     figure
%     subplot(1,2,1)
%     imagesc(abs(iq));set(gca,'YDir','normal'); colorbar;
%     subplot(1,2,2)
%     imagesc(abs(iq1));set(gca,'YDir','normal'); colorbar;
    %-------------------- POWER SPECTRUM ESTIMATION (PSE) --------------------
    siq=iq(1:B(1)*floor(size(iq,1)/B(1)),:);                             % IQ segment for PSE
    [~,h0]=iq2ps(siq,B,L);                      % zero-phase PSF

    %--------------- ZERO-PHASE VS PHASE-ADAPTIVE DECONVOLUTION --------------
    %f1=hybid(iq,h0,[lmd a],'TOT',1);            % "zero-phase" deconvolution
    [f2,h]=hybid(iq,h0,[lmd a],'TOT',100);       % full blind deconvolution

    %------------------ PSF RECONSTRUCTION VIA "MODULATION" ------------------
    psf=blur2psf(h,fn,D,C);                       % conversion to a real PSF
    %--------------Conversion to circulant matrix
    [m_hh,n_hh] = size(psf); 
    shift_h = zeros(Nz*Nx, Nt);
    shift_h(1:m_hh,1:n_hh) = psf;
    H = fft2(circshift(shift_h, 1-[floor((m_hh+1)/2),floor((n_hh+1)/2)])); 
    %-------------------------------- RESULTS --------------------------------
%     DR=20;                                      % dynamic range (DR) parameter
% 
%     env2img=@(z)uint8(round((255/log(DR+1))*... % DR normalization function
%         log((DR/max(z(:)))*z+1)));
% 
%     I0=env2img(abs(iq));                        % original envelope
%     I1=env2img(abs(f1));                        % zero-phase reconstruction
%     I2=env2img(abs(f2));                        % phase-adaptive reconstruction
% 
%     figure
%     subplot(131), imshow(I0,'Init','fit'), title('ORIGINAL')
%     subplot(132), imshow(I1,'Init','fit'), title('ZERO-PHASE')
%     subplot(133), imshow(I2,'Init','fit'), title('PHASE-ADAPTIVE')
%     set(gcf,'color','w')
% 
%     figure
%     imagesc(psf)
%     colormap(gray)
end