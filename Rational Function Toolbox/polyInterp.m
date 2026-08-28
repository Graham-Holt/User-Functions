function p = polyInterp(x,f,d)
% polyInterp(x,f) generates the simplest polynomial function defined by its 
% derivatives at various points
% 
% Graham Holt, April 2026. Updated August 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% polyInterp(x,f)
% polyInterp(___,d)
% p = polyInterp(___)
% [p,V] = polyInterp(___)
% 
%% Description
% polyInterp(x,f) returns the simplest polynomial which evaluates to "f" at
% "x"
%
% polyInterp(___,d) considers the order "d" derivative of the polynomial to
% equal "f" at "x".

% Ensures that there exists an interpolating solution (no least-squares)
if nargin<3
    d = zeros(length(x),1);
end
if ~(length(d)==length(x) && length(f)==length(x))
    error('Vetcor inputs must have the same length');
end
if any(d>=length(d)) || any(d~=round(d))
    error('Order of derivative must be an integer less than or equal to the order of polynomial.');
end
x = reshape(x,[],1); f = reshape(f,[],1);

% Initializes Vandermonde matrix
n = length(x); V = zeros(n);
if ~(isa(x,'numeric')&isa(f,'numeric'))
    V = sym(V);
end

% Constructs matrix and augments with function values
V = cvander(x,[],d); Vaug = rref([V f]); r = rank(V);

if any(Vaug((r+1):n,n+1)~=0)
    error(['(',num2str(sum(Vaug((r+1):n,n+1)~=0)),') conditions are incompatible.']);
end

% Finds minimum norm solution
V = Vaug(1:r,1:n);
f = Vaug(1:r,n+1);
p = V.'/(V*V.')*f;

% Reduces polynomial to minimal degree
p = p(find(p~=0,1):end).';

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