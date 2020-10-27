function [F,G] = sbasis(N,M)

% SBASIS Periodic spline basis
%   SBASIS(N,M) computes a discrete shift-invariant basis of M cubic
%   b-splines of length N. Note that the splines are mutually shifted
%   by interger multiples of T=N/M, which is therefore expected to be
%   an integer.
%                           [F,G] = sbasis(N,M)
% Input:
%                   N - length of splines
%                   M - number of splines (M>=4)
% Output:
%                   F - spline basis
%                   G - dual (spline) basis (i.e., G'*F = F'*G = I)
%
% written by Oleg Michailovich, 2006
% revised by Oleg Michailovich, April, 2014

if mod(N,M)~=0
    error('Unacceptable dimensions!');
else
    T=N/M;
end

n=(0:N-1)'/T;

i1=(n<=0);
i2=((0<n)&(n<=1));
i3=((1<n)&(n<=2));
i4=((2<n)&(n<=3));
i5=((3<n)&(n<=4));
i6=(n>4);

f1=0;
f2=(1/6)*(n.^3);
f3=(1/6)*(n.^3)-(4/6)*(n-1).^3;
f4=(1/6)*(-n+4).^3-(4/6)*(-n+3).^3;
f5=(1/6)*(-n+4).^3;
f6=0;

f=f1.*i1+f2.*i2+f3.*i3+f4.*i4+f5.*i5+f6.*i6;
f=circshift(f,-2*T);

F=zeros(N,M);
for k=1:M
    F(:,k)=circshift(f,(k-1)*T);
end

if (nargout>1)
    G=F/(F'*F);
end

end