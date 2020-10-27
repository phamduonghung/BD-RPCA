function [H,h] = filter2D(n,varargin)

% FILTER2D Frequency domain filter design
%   H = FILTER2D(n,c,...) computes the frequency response of a 2-D low-
%   pass filter of size 'n', with normalized cut-off (-3dB) frequencies
%   defined by 'c' (with 'c = 1' correspoinding to the full band).
%
%                       [H,h] = filter2D(n,c,e,l)
%   Input:
%                   n - size of the filter
%                   c - cut-off frequencies (default:  c = 0.5)
%                   e - roll-off rate (default: e = 0.001)
%                   l - half-size of impulse response
%
%   FILTER2D(...,'type',TYPE), with 'TYPE' being either 'sep' or 'nsep',
%   allows one to choose between separable & nonseparable filter design,
%   respectively.
%
%   [H,h] = FILTER2D(...) also computes the impulse response of the fil-
%   ter using windowing desing with a Hamming window.
%
%   Output:
%                   H - frequency response
%                   h - impulse response
%
%   See also GKERNEL, FILTS
%
% written by Oleg Michailovich, 2014/03
% revised by Oleg Michailovich, 2017/06/20
% revised by Oleg Michailovich, 2018/07/06

val4n=@(x)validateattributes(x,{'double'},{'integer'});
val4c=@(x)validateattributes(x,{'double'},{'>',0,'<',1});
val4e=@(x)validateattributes(x,{'double'},{'>',0});
val4flag=@(x)any(validatestring(x,{'sep','nsep'}));

p=inputParser;
addRequired(p,'n',val4n);
addOptional(p,'c',0.5,val4c);
addOptional(p,'e',1e-3,val4e);
addOptional(p,'l',0,val4n);
addParameter(p,'type','nsep',val4flag);
parse(p,n,varargin{:});

n=p.Results.n;
if isscalar(n)
    n=n([1 1]);
end

c=p.Results.c;
if isscalar(c)
    c=c([1 1]);
end

e=p.Results.e;

c=pi*c;
an=log(2+1/e)/(2*c(1)^2);
am=log(2+1/e)/(2*c(2)^2);
wn=(2*pi/n(1))*(-floor(n(1)/2):ceil(n(1)/2)-1);
wm=(2*pi/n(2))*(-floor(n(2)/2):ceil(n(2)/2)-1);
Hn=exp(-an.*(wn.*wn));
Hm=exp(-am.*(wm.*wm));

switch p.Results.type
    case 'nsep'
        H=(Hn'*Hm).^2;
        H=(1+e)*(H./(H+e));
    case 'sep'
        Hn2=Hn.*Hn;
        Hn2=Hn2./(Hn2+e);        
        Hm2=Hm.*Hm;
        Hm2=Hm2./(Hm2+e);        
        H=(1+e)^2*(Hn2'*Hm2);
end

H=ifftshift(H);

if (nargout>1)
    if isequal(p.Results.l,0)
        h=real(ifftshift(ifft2(H)));
    else
        if isscalar(p.Results.l)
            l=p.Results.l([1 1]);
        else
            l=p.Results.l;
        end
        
        m=2*l+1;        
        h=circshift(ifft2(H),l);
        h=win2d(m).*real(h(1:m(1),1:m(2)));
    end
end

end