function P = polyInterp(x,f,d,tol)
% polyInterp(x,f) generates the simplest polynomial function defined by its 
% derivatives at various points
% 
% Graham Holt, April 2026. Updated September 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% polyInterp(x,f)
% polyInterp(___,d)
% P = polyInterp(___)
% 
%% Description
% polyInterp(x,f) returns the simplest polynomial which evaluates to "f" at
% "x"
%
% polyInterp(___,d) considers the order "d" derivative of the polynomial to
% equal "f" at "x".
%
% polyInterp(___,tol) allows you to adjust the tolerance for the polynomial
% to be considered an interpolant (default: 1e-6)

% Ensures that there exists an interpolating solution (no least-squares)
if nargin < 3
    d = zeros(length(x),1);
end
if nargin < 4
    tol = 1e-6;
end
if ~(length(d)==length(x) && length(f)==length(x))
    error('Vetcor inputs must have the same length');
end
if any(d>=length(d)) || any(d~=round(d))
    error('Order of derivative must be an integer less than or equal to the order of polynomial.');
end
x = reshape(x,[],1); 
f = reshape(f,[],1);

% Finds minimum norm solution
V = cvander(x,[],d);
if rank(V) ~= rank([V, f])
    error(['(',num2str(rank([V, f]) - rank(V)),') conditions are incompatible.']);
end

P = f.'/(V*V.')*V;

% Reduces polynomial to minimal degree
P = polyShorten(P);
if any(abs(polyvalDiff(P,x,d)-f)>tol)
    error('Conditions did not converge');
end

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

function P = polyvalDiff(p,x,d)
% polvalDiff(p,x) evaluates polynomials and their derivatives.
% 
% Graham Holt, October 2026. Updated October 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% polvalDiff(p,x)
% polvalDiff(p,x,d)
% P = polvalDiff(___)
% 
%% Description
% polvalDiff(p,x) returns the same result as the usual polyval() function
%
% polvalDiff(p,x,d) returns a similar result with the value of the (d)th
% derivative of the polynomial at (x)

if nargin < 3
    d = zeros(1,length(x));
end
if length(d) ~= length(x)
    error('(d) must have the same length as (x)');
end

for k = 1:length(x)
    P(k,:) = polyval(polyDiffInt(p,d(k)),x(k));
end

end