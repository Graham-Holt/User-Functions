function [N,D] = ratSimplify(N,D,tol)
% ratSimplify(n,d) eliminates common factors in the numerator and 
% denominator of a rational function
% 
% Graham Holt, August 2026. Updated August 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% ratSimplify(n,d)
% ratSimplify(___,tol)
% [N,D] = ratSimplify(___)
% 
%% Description
% ratSimplify(n,d) returns two row vectors with the coefficients of the 
% simplified rational function

if ~exists('tol','var')
    tol = 1e-6;
end

rootD = polyRoots(D,tol);

for k = 1:(length(D)-1)
    n = length(N);

    N = polyReduc(N,rootD(k));
    if length(N)~=n
        D = polyReduc(D,rootD(k));
    end
end