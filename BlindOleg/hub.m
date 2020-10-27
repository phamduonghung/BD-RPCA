function [x] = hub(x,a)

% HUB Huber function
%   HUB(x,a) returns the values of the Huber function (aka the Moreau
%   envelope of |x|) at 'x', computed according to
%
%                       / |x|^2,        if |x| <= a
%               h(x) = <
%                       \ 2*a*|x|-a^2,  otherwise.
%
%   See also PROXHUB, PROXES
%
% written by Oleg Michailovich, 2014/04/29
% revised by Oleg Michailovich, 2017/05/24

if isscalar(a)
    a=a(ones(size(x)));
end

xa=abs(x);
I=(xa<=a);
J=~I;

x(I)=xa(I).^2;
x(J)=2*(xa(J).*a(J))-a(J).^2;

end