function [f,e] = proxgrad(P,NIT)

% PROXGRAD Proximal gradient method
%   PROXGRAD(f,...) implements the Proximal Gradient Method (PGM) with 
%   Nesterov acceleration (aka FISTA).
%
%                       [f,e] = proxgrad(P,NIT)
%   Input:
%                   P - structure containing:
%                           * P.f   - initialization
%                           * P.cst - cost function
%                           * P.grd - (scaled) gradient of the differen-
%                                     tiable part of P.cst
%                           * P.prx - (scaled) proximal operator of the
%                                     non-differentiable part of P.cst
%                 NIT - total number of iterations
%   Output:
%                   f - solution (either a local or global minimizer)
%                   e - cost function values
%
%   See also PROXES
%
% Written by Oleg Michailovich, 2018/07/17

[fold,fhat]=deal(P.f);
told=1;

if (nargout>1)
    e=zeros(NIT,1);
    for itr=1:NIT
        % fprintf('iteration %d of %d\n',[itr NIT]);
        
        f=P.prx(fhat-P.grd(fhat));
        e(itr)=P.cst(f);
        
        if (itr>1)&&(e(itr)>e(itr-1))
            told=1;
            fhat=fold;
            continue
        end
        
        t=(1+sqrt(1+4*told^2))/2;
        fhat=f+((told-1)/t)*(f-fold);
        fold=f;
        told=t;
    end
else
    for itr=1:NIT
        f=P.prx(fhat-P.grd(fhat));
        
        t=(1+sqrt(1+4*told^2))/2;
        fhat=f+((told-1)/t)*(f-fold);
        fold=f;
        told=t;
    end
end

end