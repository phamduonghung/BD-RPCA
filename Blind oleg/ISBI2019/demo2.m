clc
clear
close all

%------------------------- ESTIMATION PARAMETERS -------------------------
C=[0.03 0.2];                               % LPF cut-off frequencies
B=[24 16];                                  % block size for PSE
M=[8 8];                                    % #splines for PSE

lmd1=2e-3;                                  % reg. parameter (1st harm.)
lmd2=4e-3;                                  % reg. parameter (2nd harm.)
a=0.05;                                     % Huber parameter

%-------------------- DATA READING AND NORMALIZATION  --------------------
d=load('data.mat');                         % reading RF data

%----------------------------- DEMODULATION ------------------------------
fn1=d.f0./d.fs;                             % norm. cent. freq. (1st)
fn2=2*fn1;                                  % norm. cent. freq. (2nd)

iq1=rf2iq(d.rf,fn1,C);                      % demodulation (1st harm.)
iq2=rf2iq(d.rf,fn2,C);                      % demodulation (2nd harm.)

iq1=iq1/(2^ceil(log2(max(abs(iq1(:))))));   % normalization
iq2=iq2/(2^ceil(log2(max(abs(iq2(:))))));   % normalization

%-------------------- POWER SPECTRUM ESTIMATION (PSE) --------------------
[~,h01]=iq2ps(iq1(1:96,:),B,M);             % initial PSF (1st harm.)
[~,h02]=iq2ps(iq2(1:96,:),B,M);             % initial PSF (2nd harm.)

[f1,h1]=hybid(iq1,h01,[lmd1 a],'TOT',20);   % blind deconv. (1st harm.)
[f2,h2]=hybid(iq2,h02,[lmd2 a],'TOT',20);   % blind deconv. (2nd harm.)
 
%-------------------------------- RESULTS --------------------------------
DR=20;                                      % dynamic range (DR) parameter

env2img=@(z)uint8(round((255/log(DR+1))*... % DR normalization function
    log((DR/max(z(:)))*z+1)));

I01=env2img(abs(iq1));                      % envelope (1st harmonic)
I02=env2img(abs(iq2));                      % envelope (2nd harmonic)
I1=env2img(abs(f1));                        % 1st harmonic reconstruction
I2=env2img(abs(f2));                        % 2nd harmonic reconstruction

figure
subplot(221), imshow(I01,'Init','fit'), title('ORIGINAL (1ST)')
subplot(222), imshow(I02,'Init','fit'), title('ORIGINAL (2ND)')
subplot(223), imshow(I1,'Init','fit'), title('ESTIMATED (1ST)')
subplot(224), imshow(I2,'Init','fit'), title('ESTIMATED (2ND)')
set(gcf,'color','w')