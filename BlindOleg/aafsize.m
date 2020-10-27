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