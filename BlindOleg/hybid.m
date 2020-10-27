function [trf,psf] = hybid(g,h,varargin)

% HYBID "Hybrid" blind image deconvolution of IQ images.
%   HYBID(g,h,...) deconvolves an IQ image 'g' using a "zero-phase" version
%   'h' of the true (yet unknown) point-spread function (PSF). Specifically,
%   the deconvolution solves the following optimization problem:
%
%               min_{f,s} {||f * s - g||^2 + lmd * phi(f)},
%                           s.t. abs(fft2(s)) = H,
%
%   with 'f' standing for the tissue reflectivity function to be recovered,
%   H=abs(fft2(h)) being the magnitude spectrum of the PSF, and 'phi' being
%   a regularization functional.
%
%   HYBID(g,h,prm) uses 'prm' to control the type of functional 'phi' along
%   with the value of its associated regularization parameter 'lmd'. In par-
%   ticular, while the latter is always given by 'prm(1)', the type of 'phi'
%   depends on whether 'prm' is a scalar or a length-2 vector. In the former
%   case, 'phi' is set to be the l1-norm, whereas in the latter case, it is
%   defined as Huber function with a smoothing parameter 'a = prm(2)'.
%
%   HYBID(...,'PropertyName',PropertyValue,...) can be used to control some
%   optimization parameters, as described in the table below.
%
%      *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*
%      | PropertyName | PropertyValue |          Description          |
%      |              |   [default]   |                               |
%      *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*
%      |      TOT     |       10      | Maximum #interations of AM.   |          |
%      |      NIT     |       50      | Maximum #interations of PGM.  |
%      |      tol     |      5e-3     | Tolerance for AM termination. |
%      |      vrb     |      True     | Verbosity (on/off)            |
%      |      all     |      False    | Output format (see below)     |
%      *--------------------------------------------------------------*
%
%   [trf,psf] = HYBID(...) returns the estimates of tissue reflectivity and
%   PSF. When 'all' is set to False (default), only the final solutions are
%   returned. Alternatively, when 'all' is set to True, the function returns
%   all intermediate estimates in the forms of cell arrays.
%
%   See also KILLEDGE, PROXES, HUB, PROXGRAD
%
%   Written by Oleg Michailovich, 2018/10/23.

ifmtx=@(x)validateattributes(x,{'double'},{'2d'});
ifvec=@(x)validateattributes(x,{'numeric'},{'vector','nonnegative'});
ifint=@(x)validateattributes(x,{'numeric'},{'scalar','positive','integer'});
ifrea=@(x)validateattributes(x,{'numeric'},{'scalar','positive','real'});
iflog=@(x)validateattributes(x,{'numeric'},{'scalar','binary'});

p=inputParser;
p.CaseSensitive=true;
addRequired(p,'g',ifmtx);
addRequired(p,'h',ifmtx);
addOptional(p,'lmd',[0.1 0.05],ifvec);
addParameter(p,'TOT',10,ifint);
addParameter(p,'NIT',50,ifint);
addParameter(p,'tol',5e-3,ifrea);
addParameter(p,'all',false,iflog);
addParameter(p,'vrb',true,iflog);
parse(p,g,h,varargin{:});

N=size(p.Results.h);
L=size(p.Results.g);

h=p.Results.h;
W=abs(fft2(h.*win2d(N),L(1),L(2)));
W=(W/max(W(:))).^2;
W=sqrt(W./(W+1e-2));

P.f=zeros(L+N-1);
g=killedge(p.Results.g,h);
G=fft2(g);

if isscalar(p.Results.lmd)
    P.prx=@(x)(proxes('s',x,p.Results.lmd));
    RG=@(x)sum(abs(x(:)));
    reg='l1-norm';
else
    P.prx=@(x)(proxes('h',x,p.Results.lmd(1),p.Results.lmd(2)));
    RG=@(x)sum(hub(x(:),p.Results.lmd(2)));
    reg=['Huber (a=',mat2str(p.Results.lmd(2)),')'];
end

if p.Results.all
    trf=cell(p.Results.TOT,1);
    psf=cell(p.Results.TOT,1);
end

if p.Results.vrb
    fprintf('\nAM iterations are initialized with:\n')
    fprintf('Max number of AM iterations: \t%d\n',p.Results.TOT)
    fprintf('Number of PGM iterations: \t%d\n',p.Results.NIT)
    fprintf('Termination tolerance: \t\t%4.2f\n',p.Results.tol)
    fprintf('Regularization parameter: \t%6.4f\n',p.Results.lmd(1))
    fprintf('Regularization method: \t\t%s\n',reg)
    fprintf('\n')
end

vct=@(z)z(:);

for itr=1:p.Results.TOT
    hh=conj(flip(flip(h,1),2));
    HH=@(x)conv2(x,h,'valid');
    HA=@(x)conv2(x,hh);
    ER=@(x)(HH(x)-g);
    
    P.cst=@(x)(0.5*sum(vct(abs(ER(x)).^2))+p.Results.lmd(1)*RG(x));
    P.grd=@(x)HA(ER(x));
    
    P.f=proxgrad(P,p.Results.NIT);
    cst=P.cst(P.f);
    if p.Results.vrb
        fprintf('AM iter %d out of %d, cost: %f\n',[itr p.Results.TOT cst])
    end
    
    h=phasadj(P.f,G,h,W);
    
    if p.Results.all
        trf{itr}=wkeep1(P.f,L);
        psf{itr}=h;
    else
        trf=wkeep1(P.f,L);
        psf=h;
    end
    
    if (itr==1)
        cst_old=cst;
    else
        cst_new=cst;
        nerr=(cst_old-cst_new)/cst_old;
        
        if (nerr>=0)&&(nerr<p.Results.tol)
            if p.Results.all
                trf=trf(1:itr);
                psf=psf(1:itr);
            end
            
            break
        else
            cst_old=cst_new;
        end
    end
end

end

function [h] = phasadj(f,G,h,W)

L=size(G);
M=numel(h);
N=size(h)-1;

if nargin<4
    W=1;
end

w1=(2*pi/L(1))*(0:L(1)-1);
w2=(2*pi/L(2))*(0:L(2)-1);
[w1,w2]=ndgrid(w1,w2);
[n1,n2]=ndgrid(0:N(1),0:N(2));
Phi0=w1(:)*n1(:)'+w2(:)*n2(:)';

P=G.*conj(fft2(conv2(f,h,'valid')));
PhiR=-N(1)*w1-N(2)*w2-angle(P);
PhiR=atan2(sin(PhiR),cos(PhiR));
PhiR=(PhiR+(2*pi)*(PhiR<0))/2;

Phi=Phi0+PhiR(:);
A=[sin(Phi) -cos(Phi)].*W(:);
[r,~]=eigs(A'*A,1,'SM');
r=reshape(complex(r(1:M),r(M+1:end)),N+1);

denom=fft2(r);
num=fft2(conj(flip(flip(r,1),2)));
h=ifft2(fft2(h).*exp(1i*angle(num./denom)));

end