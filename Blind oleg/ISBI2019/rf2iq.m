function [iq,D] = rf2iq(rf,fn,C,D)

% RF2IQ Demodulation of ultrasound RF images.
%   RF2IQ(rf,...) converts an RF image 'rf' into its associated IQ-image using
%   frequency-domain demodulation and anti-aliasing filtering.
%
%                         [iq] = rf2iq(rf,fn,C,D)
%   Input: 
%              rf - RF image
%              fn - normalized central frequency
%               C - normalized cut-off frequencies
%               D - decimation rates (if omitted, D is computed from C) 
%   Output:
%              iq - IQ image
%               D - decimation rates
%
%   See also FILTER2D
%
%   Written by Oleg Michailovich, 07/2018.
%   Revised by Oleg Michailovich, 10/2018.
%rf=X;
if ~ismatrix(rf)
    error('The first argument must be a matrix')
end

if isscalar(C)
    C(2)=1-eps;
end

n=size(rf);
N=aafsize(filter2D(n,C,1e-3));
N=N-rem(N,2)+1;
M=(N-1)/2;

L=n+N-1;
T=-round(L(1)*fn);
B=filter2D(L,C,1e-3);

F=fft2(rf,L(1),L(2));
iq=circshift(ifft2(B.*circshift(F,T)),M);
iq=iq(M(1)+1:M(1)+n(1),M(2)+1:M(2)+n(2));

if 0
    figure
    subplot(2,3,1)
    imagesc(abs(fftshift(B))); set(gca,'YDir','normal'); colorbar; title('filter B')
    % 
    % F=fft2(rf,L(1),L(2));
    subplot(2,3,2)
    imagesc(abs(fftshift(F))); set(gca,'YDir','normal'); colorbar; title('fft(rf)')
    F1=fft2(real(rf),L(1),L(2));
    subplot(2,3,3)
    imagesc(abs(fftshift(F1))); set(gca,'YDir','normal'); colorbar; title('fft real(rf)')
    % 
    % iq=circshift(ifft2(B.*circshift(F,T)),M);
    subplot(2,3,4)
    imagesc(abs(fftshift(B.*circshift(F,T)))); set(gca,'YDir','normal'); colorbar; title('filter B*filter')
    subplot(2,3,5)
    imagesc(abs(iq)); set(gca,'YDir','normal'); colorbar; title('iq')
    % 
    % iq=iq(M(1)+1:M(1)+n(1),M(2)+1:M(2)+n(2));
    subplot(2,3,6)
    imagesc(abs(iq)); set(gca,'YDir','normal'); colorbar; title('undersampling iq')
end

if nargin<5
    D=max(floor(1./C-1),1);
    iq=iq(1:D(1):n(1),1:D(2):n(2));
else
    if iscalar(D)
        D(2)=1;
    end
    iq=iq(1:D(1):n(1),1:D(2):n(2));
end

end

function [NAF] = aafsize(B)

n=size(B);
x=(0.5-n(1)/2:n(1)/2-0.5)';
y=(0.5-n(2)/2:n(2)/2-0.5)';

b=abs(fftshift(ifft2(B)));

b1=sum(b,2); 
b1=b1(:)/sum(b1);
b2=sum(b,1);
b2=b2(:)/sum(b2);

s1=sqrt(b1'*(x-b1'*x).^2);
s2=sqrt(b2'*(y-b2'*y).^2);
NAF=floor(2.5*[s1 s2]);

end