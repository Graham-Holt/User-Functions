function [N,D] = ratInterp(x,f,asym)
% ratInterp(x,f,asym) generates the simplest interpolating rational 
% function with given asymptotic degree
% 
% Graham Holt, September 2026. Updated September 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% ratInterp(x,f)
% ratInterp(___,asym)
% p = ratInterp(___)
% 
%% Description
% polyInterp(x,f) returns the simplest polynomial which evaluates to "f" at
% "x"
%
% polyInterp(___,asym) considers the degree of the asymptote for the
% rational function.

% Ensures that there exists an interpolating solution (no least-squares)
if nargin<3
    asym = length(x) - 1;
end
if ~(length(f)==length(x))
    error('Vetcor inputs must have the same length');
end
if asym >= length(x) || asym~=round(asym) || mod(asym,2)~=mod(length(x)-1,2)
    error('Degree of asymptote must be an integer less than or equal to the order of rational and of the same parity.');
end
x = reshape(x,[],1); 
f = reshape(f,[],1);
n = 0.5*(length(x) + asym - 1);
m = 0.5*(length(x) - asym - 1);

% Constructs matrix and augments with function values
V = [cvander(x,n+1), -f.*cvander(x,m+1)];

% Finds minimum norm solution and simplifies coefficients
coeff = null(V).'; I = find(abs(coeff)>1e-12,1,'last'); 

if all(coeff((end-m):end) <= 1e-12)
    error('Conditions are incompatible.');
end

N = coeff(1:(n+1))/coeff(I);
N = N(find(abs(N)>1e-12,1):end);

D = coeff((end-m):end)/coeff(I);
D = D(find(abs(D)>1e-12,1):end);

end

function V = cvander(v,nCols,d)
% cvander(v) computes confluent Vandermonde matrices.
% 
% Graham Holt, August 2026. Updated August 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% cvander(v)
% cvander(v,nCols)
% V = cvander(___)
% 
%% Description
% cvander(v) returns a square Vandermonde matrix with the values in "v"
%
% cvander(v,nCols) returns a rectangular matrix with "nCol" columns
% 
% cvander(___,d) applies the given order of derivative to the corresponding
% row of the rectangular Vandermonde matrix.

if nargin < 2
    nCols = length(v);
end
if isempty(nCols)
    nCols = length(v);
end
if nargin < 3
    d = zeros(1,length(v));
end

del = nCols - (1:nCols);
for k = 1:length(v)
    V(k,:) = [v(k).^del((1+d(k)):nCols) zeros(1,d(k))];

    for j = 1:d(k)
        V(k,:) = V(k,:).*circshift(del,1-j);
    end
end

end