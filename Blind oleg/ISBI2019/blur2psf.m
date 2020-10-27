function [s] = blur2psf(h,fn,D,C)

% BLUR2PSF Complex-to-real PSF conversion.
%   BLUR2PSF(h,...) converts a complex-valued point-spread function (PSF) 'h' 
%   into its real-valued conterpart.
%
%                       [s] = blur2psf(hc,fn,D)
%   Input:
%                  hc - complex PSF 
%                  fn - normalized central frequency
%                   D - axial [D(1)] and lateral [D(2)] decimation rates.
%                       If D is a scalar, then D(2) is set to 1.
%   Output: 
%                   s - real-valued PSF
%
% Written by O. Michailovich, 2018/10/26.

if isscalar(D)
    D=[D 1];
end

P=floor((D-1)/2);
hh=upsample(upsample(h,D(1),P(1)).',D(2),P(2)).';

n=size(hh); 
m=2.^ceil(log2(1.5*n));
H=filter2D(m,C,1e-2);

HH=zeros(m); 
HH(1:n(1),1:n(2))=hh;
HH=circshift(fft2(HH).*H,[round(fn*m(1)) 0]);
HH(2:m(1),2:m(2))=HH(2:m(1),2:m(2))+conj(HH(m(1):-1:2,m(2):-1:2));
HH(1,2:m(2))=HH(1,2:m(2))+conj(HH(1,m(2):-1:2));
HH(2:m(1),1)=HH(2:m(1),1)+conj(HH(m(1):-1:2,1));
HH(1)=0;

s=ifft2(HH./max(abs(HH(:))));
s=s(1:n(1),1:n(2));

end