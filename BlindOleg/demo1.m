clc
clear
close all

%------------------------- ESTIMATION PARAMETERS -------------------------
C=[0.03 0.2];                               % LPF cut-off frequencies
B=[24 16];                                  % block size for PSE
M=[8 8];                                    % #splines for PSE

lmd=2e-3;                                   % regularization parameter
a=0.05;                                     % Huber parameter

%-------------------- DATA READING AND NORMALIZATION  --------------------
d=load('data.mat');                         % reading RF data
[~,d.rf]=dynarange(d.rf);                   % normalization of RF data

%----------------------------- DEMODULATION ------------------------------
fn=d.f0./d.fs;                              % normalized central frequency
[iq,D]=rf2iq(d.rf,fn,C);                    % demodulation of RF data

%-------------------- POWER SPECTRUM ESTIMATION (PSE) --------------------
siq=iq(1:96,:);                             % IQ segment for PSE
[~,h0]=iq2ps(siq,B,M);                      % zero-phase PSF

%--------------- ZERO-PHASE VS PHASE-ADAPTIVE DECONVOLUTION --------------
f1=hybid(iq,h0,[lmd a],'TOT',1);            % "zero-phase" deconvolution
[f2,h]=hybid(iq,h0,[lmd a],'TOT',10);       % full blind deconvolution

%------------------ PSF RECONSTRUCTION VIA "MODULATION" ------------------
psf=blur2psf(h,fn,D);                       % conversion to a real PSF

%-------------------------------- RESULTS --------------------------------
DR=20;                                      % dynamic range (DR) parameter

env2img=@(z)uint8(round((255/log(DR+1))*... % DR normalization function
    log((DR/max(z(:)))*z+1)));

I0=env2img(abs(iq));                        % original envelope
I1=env2img(abs(f1));                        % zero-phase reconstruction
I2=env2img(abs(f2));                        % phase-adaptive reconstruction

figure
subplot(131), imshow(I0,'Init','fit'), title('ORIGINAL')
subplot(132), imshow(I1,'Init','fit'), title('ZERO-PHASE')
subplot(133), imshow(I2,'Init','fit'), title('PHASE-ADAPTIVE')
set(gcf,'color','w')

figure
imagesc(psf)
colormap(gray)