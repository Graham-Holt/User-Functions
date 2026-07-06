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
% polyMod(___,tol) uses the value "tol" to compute and distinguish roots to 
% this level of precision. It is recommended to use the default (1e-6) if 
% you do not have Symbolic Math Toolbox.

if ~exist('tol','var')
    tol = 1e-6;
end

% Evaluates roots to choosen level of detail
n = length(m)-1;

if license('test','symbolic_toolbox')
    syms x; d0 = digits; digits(ceil(-log10(tol)));
    r = vpasolve((x.^flip(0:n))*(m.')==0,x); r = double(r);
    digits(d0);
else
    r = round(roots(m),ceil(-log10(tol)));
end
% Determines derivative needed to accomodate repeated roots
k = 1; d = zeros(1,n);
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

% Compute reduced polynomial
p = polyInterp(r,d,f);

end