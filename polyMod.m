function p = polyMod(p,m,tol)
% polyMod(p,m) finds a reduction of the given polynomial with respect to a
% polynomial modulus.
% 
% Graham Holt, July 2026. Updated July 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% polyMod(p,m)
% polyMod(___,tol)
% p = polyMod(___)
% 
%% Description
% polyMod(p,m) returns the reduced polynomial of "p" modulo "m".
%
% polyMod(___,tol) uses the value "tol" to distinguish roots to this level
% of precision. It is recommended to use the default 1e-6, as this is the
% general precision of the roots() function

if ~exist('tol','var')
    tol = 1e-6;
end

% Evaluates roots to choosen level of detail
r = round(roots(m),ceil(-log10(tol)));

% Determines derivative needed to accomodate repeated roots
n = length(r); k = 1; d = zeros(1,n);
while k<=n
    t = sum(abs(r-r(k))<=tol);
    d(k:(k+t-1)) = 0:(t-1);
    k = k + t;
end

% Evaluate polynomial derivatives at modulus roots
del = length(p) - (1:length(p));
for k = 1:n
    pTemp = p;
    for j = 1:d(k)
        pTemp = circshift(pTemp.*del,1);
    end
    f(k) = polyval(pTemp,r(k));
end

end

% Compute reduced polynomial
p = polyInterp(r,d,f);

function p = polyInterp(x,d,f)
% polyInterp(x,d,f) generates the simplest polynomial function defined by 
% its derivatives at various points
% 
% Graham Holt, April 2026. Updated July 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% polyInterp(x,d,f)
% p = polyInterp(___)
% 
%% Description
% polyInterp(x,d,f) returns the simplest polynomial which has the 
% derivatives "d" equal to "f" at "x"

% Ensures that the polynomial is possible
if ~(length(d)==length(x) && length(f)==length(x))
    error('Vetcor inputs must have the same length');
end
if any(d >= length(d)) || any(d~=round(d))
    error('Order of derivative must be an integer less than or equal to the order of polynomial.');
end
x = reshape(x,[],1); f = reshape(f,[],1);

% Initializes matrix for solution
n = length(x); P = zeros(n);
if ~(isa(x,'numeric')&isa(f,'numeric'))
    P = sym(P);
end

% Constructs matrix and augments with function values
del = n - (1:n);
for k = 1:n
    P(k,:) = x(k).^del;
    for j = 1:d(k)
       P(k,:) = circshift(del.*P(k,:),1);
    end
end
Paug = rref([P f]); r = rank(P);

% Finds singular or minimum norm solution
if r==n 
    p = P\f;
elseif r~=n && Paug(end,end)==0
    P = Paug(1:r,1:(end-1));
    f = Paug(1:r,end);
    p = P.'/(P*P.')*f;
else
    error('One or more conditions are incompatible.');
end

% Reduces polynomial to minimal degree
p = p(find(p~=0,1):end).';
end