function [H,h] = iq2ps(iq,B,L,symm)

% IQ2PS Estimation of the Fourier magnitude of PSF.
%   IQ2PS(iq,B,...) estimates the Fourier magnitude of the point-spread
%   function (PSF) associated with IQ-image 'iq'. The estimation starts
%   with dividing 'iq' into blocks of size 'B(1)xB(2)' and then follows
%   the homomorphic denoising method of O. Michailovich and D. Adam.
%
%   Input:
%                  iq - IQ-image
%                   B - size of IQ-blocks
%                   L - number of b-splines along each coordinate
%                symm - if true, lateral symmetry is enforced (default)
%
%   Output:
%                   H - estimated Fourier magnitude
%                   h - zero-phase PSF
%
%   See also SBASIS, IM2BLOCKS, WIN2D
%
%   Written by O. Michailovich, 07/2018.
%   Revised by O. Michailovich, 20/10/2018.

if nargin<4
    symm=true;
end

if isscalar(B)
    B=[B B];
end

if isscalar(L)
    L=[L L];
end

B=ceil(B./L).*L;
Z=im2blocks(iq,B);
[U1,V1]=sbasis(B(1),L(1));
[U2,V2]=sbasis(B(2),L(2));

win=win2d(B);
s2n=@(x)(mean(abs(x(:)))/std(abs(x(:))));
lns=@(x)fftshift(log(abs(fft2(x.*win))));
med=@(x)fftshift(medfilt2(x,[3 3],'sym'));

w=cellfun(s2n,Z);
w=exp(-2*(w-1.91).^2);
w=permute(w,[1 3 2])/sum(w);

Z=cellfun(lns,Z,'Un',0);
Z=cellfun(med,Z,'Un',0);
H=sum(cat(3,Z{:}).*w,3);
H=exp(U1*(V1'*H*V2)*U2');
H=H/max(H(:));

if symm
    H=[H(:,1) (H(:,2:B(2))+H(:,B(2):-1:2))/2];
end

if (nargout>1)
    h=ifft2(H);
    I=(fftshift(abs(h)/h(1))>=0.02);
    J1=any(I,2);
    J2=any(I,1);
    N(1)=find(J1,1,'last')-find(J1,1,'first')+1;
    N(2)=find(J2,1,'last')-find(J2,1,'first')+1;
    N=N-rem(N,2)+1;
    
    h=circshift(h,(N-1)/2);
    h=h(1:N(1),1:N(2));
end

end