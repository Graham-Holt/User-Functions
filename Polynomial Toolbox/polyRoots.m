function R = polyRoots(p,tol)
% polyRoots(p,tol) computes polynomial roots to arbitrary precision
% 
% Graham Holt, August 2026. Updated August 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% polyRoots(p,tol)
% R = polyRoots(___)
% 
%% Description
% polyRoots(p,tol) returns a column vector with the roots of the input
% polynomial to the given degree of precision

R = roots(p);

if license('test','symbolic_toolbox')

d0 = digits; digits(ceil(-log10(tol)));
R = vpasolve(poly2sym(p)==0); R = double(R);
digits(d0);

else

dp = polyDiffInt(p,1);
res = -polyval(p,R)./polyval(dp,R);
while any(abs(res)>tol)
    R = R + res;
    res = -polyval(p,R)./polyval(dp,R);
end
R = round(R,ceil(-log10(tol)));

end