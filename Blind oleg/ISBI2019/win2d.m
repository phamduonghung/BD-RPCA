function [w] = win2d(N)

% WIN2D Two-dimensional Hanning window
%   WIN2D(N) computes a 2-D Hanning window of size N(1)-by-N(2) as an 
%   outer (aka tensor) product of two 1-D Hanning windows of respecti-
%   ve sizes.
%
%   See also WINDOW, HANNING
%
% Written by Oleg Michailovich, 2017/08/10

if isscalar(N)
    validateattributes(N,{'numeric'},...
        {'integer','positive'},'win2d','N')
    N=N([1 1]);
else
    validateattributes(N,{'numeric'},...
        {'integer','positive','size',[1 2]},'win2d','N')
end

w1=hanning(N(1));
w2=hanning(N(2));
w=w1*w2';

end