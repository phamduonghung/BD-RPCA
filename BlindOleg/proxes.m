function [y] = proxes(fun,x,lmd,a)

% PROXES Proximal mappings
%   PROXES(fun,x..) subjects the values in 'x' the proximal mapping corres-
%   ponding to the function defined by 'fun'. At this time, PROXES supports
%   the following functions.
%
%               .------------------------------------------.
%               |     fun     |      Related function      |
%               |------------------------------------------|
%               |     's'     |           lmd*|x|          |
%               |------------------------------------------|
%               |     'r'     |    (lmd/2)*||x| - a|^2     |
%               |------------------------------------------|
%               |     'l'     |     -(lmd/2)*log(x^2)      |
%               |------------------------------------------|
%               |     'h'     |        	lmd*hub(x,a)       |
%               .------------------------------------------.
%
%                       [y] = proxes(fun,x,lmd,w)
%   Input:
%                   fun - selects the type of proximal mapping
%                     x - input data to be processed
%                   lmd - regularization parameter
%                     a - additional parameter to be passed into fun = 'r' 
%                           or fun = 'h' (see the table above)
%   Output:
%                     y - processed data  
%
%   See also SOFT, HUB, PROXHUB
%
% Written by Oleg Michailovich, 2016
% Revised by Oleg Michailovich, 2018/07/26

xa=abs(x);
e=x./(xa+(xa==0));

switch fun
    case 's'
        xa=max(xa-lmd,0);
    case 'r'
        theta=1./(1+lmd);
        xa=(theta.*xa+(1-theta).*a);
    case 'l'
        xa=xa/2;
        xa=xa+sqrt(xa.^2+lmd);
    case 'h'
        lmd=2*lmd;
        xa=max(xa./(1+lmd),xa-lmd.*a);
end

y=xa.*e;

end