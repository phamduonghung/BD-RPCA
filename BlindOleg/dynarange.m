function [DR,x] = dynarange(x)

% DYNARANGE Dynamic range normalization.
%   DR = DYNARANGE(x), with 'x' being a numeric array, computes the dynamic
%   range of 'x' as the smallest dyadic integer, which bounds 'abs(x)' from 
%   above.
%
%   [~,xn] = DYNARANGE(x) returns a normalized version of 'x'.
%
%   Written by Oleg Michailovich, 10/2018.

if ~isnumeric(x)
    error('The input must be numeric!')
end

DR=max(abs(x(:)));
DR=2^(ceil(log2(DR)));

if (nargout>1)
    x=(1/DR)*x;
end

end