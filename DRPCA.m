function [T, x] = DRPCA(S, H, lambda, rho, mu, eps, mu_max, tol, max_iter)

% Deconvuluted DRPCA.
% Input
% - S is a data matrix (of the size n x m) to be decomposed
%   S can also contain NaN's for unobserved values
% - H: PFS
% - lambda - regularization parameter, default = 1/sqrt(max(n,m))
% - mu - the augmented lagrangian parameter, default = 10*lambda
% - tol - reconstruction error tolerance, default = 1e-6
% - max_iter - maximum number of iterations, default = 1000
%
% Ouput
%- x is high resolution blood
% - T is tissue

    [m, n] = size(S);
    unobserved = isnan(S);
    S(unobserved) = 0;
    normX = norm(S, 'fro');

    % default arguments
    if nargin < 2
        H = ones(m,n);
    end
    if nargin < 3
        lambda = 1./sqrt(max(m,n));
    end
    if nargin < 4
        rho = 1;
    end
    if nargin < 5
        mu = 1e-6;
    end
    if nargin < 6
        eps = 2;
    end
    if nargin < 7
        mu_max = 1e6;
    end
    if nargin < 8
        tol = 1e-6;
    end
    if nargin < 9
        max_iter = 50;
    end
    
    %% initial solution
    T = zeros(m, n);
    x = zeros(m, n);
    z = zeros(m, n);
    N = zeros(m, n); % Nu 
    W = zeros(m, n); % W
    Hx = real(ifft2(H.*fft2(x)));
    
    for iter = (1:max_iter)
        % ADMM step: update T and x, z
        T = Do(rho/mu, S - Hx + (1/mu)*N);
        z = So(lambda/mu, x + (1/mu)*W);
        x1 = real(ifft2(conj(H).*fft2(S-T))) + z + (1/mu)*(real(ifft2(conj(H).*fft2(N))) - W);
        h1 = 1./(abs(H).^2 + ones(m,n));
        x = real(ifft2(h1.*fft2(x1)));
        Hx = real(ifft2(H.*fft2(x)));
        % and augmented lagrangian multiplier
        Z1 = S - T - Hx;
        Z1(unobserved) = 0; % skip missing values
        N = N + mu*Z1; 
        Z2 = x - z;
        Z2(unobserved) = 0; % skip missing values
        W = W + mu*Z2;
        mu = min(mu*eps, mu_max);
        
        err1 = norm(Z1, 'fro') / normX;
        err2 = norm(Z2, 'fro') / normX;
        if (iter == 1) || (err1 > tol) || (err2 > tol) 
            fprintf(1, 'iter: %04d\terr1: %f\terr2: %f\trank(T): %d\tcard(S): %d\n', ...
                    iter, err1, err2, rank(T), nnz(x(~unobserved)));
        end
        if (err1 < tol) && (err2 < tol)
            break; 
        end
    end
end

function r = So(tau, S)
    % shrinkage operator
    r = sign(S) .* max(abs(S) - tau, 0);
end

function r = Do(tau, S)
    % shrinkage operator for singular values
    [U, D, V] = svd(S, 'econ');
    r = U*So(tau, D)*V';
end
