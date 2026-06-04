function T = butcherIRK(p,method,s,alpha,beta)
% butcherIRK(p) generates a Butcher tableau for an implicit Runge-Kutta 
% integrator with order "p"
% 
% Graham Holt, May 2026. Updated June 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% butcherIRK(p)
% butcherIRK(p,method)
% butcherIRK(___,s,alpha,beta)
% T = butcherIRK(___)
% 
%% Description
% butcherIRK(p) generates a stable p-order implicit Runge-Kutta Butcher 
% tableau based on the Lobatto III method.
%
% butcherIRK(p,method) generates a stable p-order implicit Runge-Kutta 
% Butcher tableau based on the method described by the input "method".
% Valid inputs include: "Gauss", "Radau1", "Radau2", "Lobatto", "Custom"
%
% butcherIRK(___,s,alpha,beta) generates a stable s-stage implicit 
% Runge-Kutta Butcher tableau of at least order "p" using the free
% variables "alpha" and "beta". These inputs are only used when the
% "Custom" method is specified.
%
% Further reading:
% Solving Ordinary Differential Equations II: Stiff and Differential - 
% Algebraic Problems by Ernst Hairer and Gerhard Wanner
% Chapter 4.5, Pages 75-96, ISBN 9783662099476 (eBook)


% Default to Lobatto III method
if ~exist('method','var')
    method = 'L';
end

% Define the stage polynomial (M) and the transformed matrix (X)
if any(strcmpi(method,{'Gauss','G'}))
    p = p + mod(p,2); s = p/2; zeta = floor(0.5*(p-1)); eta = p-zeta-1;
    
    xi = 0.5./sqrt(4*(1:(s-1)).^2 - 1);
    X = diag(-xi,1) + diag(xi,1)'; X(1,1) = 0.5;

    M = sLegendre(s);

elseif any(strcmpi(method,{'Radau1','R1'}))
    p = p + 1 - mod(p,2); s = 0.5*(p+1); zeta = floor(0.5*(p-1)); eta = p-zeta-1;

    xi = 0.5./sqrt(4*(1:(s-1)).^2 - 1);
    X = diag(-xi,1) + diag(xi,1)'; X(1,1) = 0.5;
    
    M = 0.5*(sqrt((2*s - 1)/(2*s + 1))*sLegendre(s) + [0 sLegendre(s-1)]);

elseif any(strcmpi(method,{'Radau2','R2'}))
    p = p + 1 - mod(p,2); s = 0.5*(p+1); zeta = floor(0.5*(p-1)); eta = p-zeta-1;

    xi = 0.5./sqrt(4*(1:(s-1)).^2 - 1);
    X = diag(-xi,1) + diag(xi,1)'; X(1,1) = 0.5;
    
    M = 0.5*(sqrt((2*s - 1)/(2*s + 1))*sLegendre(s) - [0 sLegendre(s-1)]);
    
elseif any(strcmpi(method,{'Radau3','R3','Lobatto','L'}))
    p = p + mod(p,2); s = 0.5*(p+2); zeta = floor(0.5*(p-1)); eta = p-zeta-1;

    xi = 0.5./sqrt(4*(1:(s-1)).^2 - 1);
    X = diag(-xi,1) + diag(xi,1)'; X(1,1) = 0.5;
    
    M = [sLegendre(s-1)./(s:-1:1) 0];

elseif strcmpi(method,{'Custom'})
    k = 2*s- p; zeta = floor(0.5*(p-1)); eta = p-zeta-1;
    xi = 0.5./sqrt(4*(1:(s-1)).^2 - 1);
    X = diag(-xi,1) + diag(xi,1)'; X(1,1) = 0.5;
    if exist('beta','var')
        [zeta,eta] = size(beta); zeta = s-zeta; eta = s-eta;
        X((zeta+1):end,(eta+1):end) = beta;
    end
    
    if p<(eta+zeta+1) || p<2*zeta+1
        error('Size of BETA too large to guarantee order conditions');
    end

    if ~exist('alpha','var')
        M = [sLegendre(s-1)./(s:-1:1) 0];
    elseif length(alpha)==k+1
        M = zeros(1,s+1);
        for j = 0:k
            M = M + [zeros(1,j) alpha(j+1)*sLegendre(s-j)];
        end
    else
        error('Size of ALPHA is inconsistent with order and stages.');
    end
    
end

% Determine stages as roots of M
c = sort(roots(M),'ascend');
if any(imag(c)~=0)|| (length(unique(c))<s)
    error('Could not produce enough distinct, real time steps');
end

% Apply order conditions to obtain weights
b = (1./(1:s))/(c.^(0:(s-1)));

% Compute W-transformation matrix
[~,W] = sLegendre(0:(s-1),c);

% Confirm order conditions and form Butcher tableau
J = W'*diag(b)*W; m = max([1,eta,zeta]);
if det(W)==0 || norm(J(1:m,1:m)-eye(m))>1e-9
    error('Condition is incompatible with requirements');
end
T = round([c W*X/W; 0 b],12);

end


function [p,L] = sLegendre(n,x)
% Computes the coefficients and values of shifted Legendre polynomials of
% order "n" at "x"

p = [ones(length(n),1), zeros(length(n),max(n))];
for k = 1:length(n)
    for j = 0:n(k)
        cp = prod(((n(k)+(1:j))./(1:j)).*(n(k)-(1:j)+1)./(1:j));
        p(k,j+1) = ((-1)^(j+n(k)))*sqrt(2*n(k)+1)*cp;
    end
    p(k,:) = flip(p(k,:));
    if exist('x','var')
        L(:,k) = polyval(p(k,:),x);
    end
end

end