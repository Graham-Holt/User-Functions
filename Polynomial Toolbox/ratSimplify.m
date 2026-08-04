function [N,D] = ratSimplify(N,D)
% ratSimplify(n,d) eliminates common factors in the numerator and 
% denominator of a rational function
% 
% Graham Holt, August 2026. Updated August 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% ratSimplify(n,d)
% [N,D] = ratSimplify(___)
% 
%% Description
% ratSimplify(n,d) returns two row vectors with the coefficients of the 
% simplified rational function

rootD = polyRoots(D,1e-12);

for k = 1:(length(D)-1)
    n = length(N);

    N = polyReduc(N,rootD(k));
    if length(N)~=n
        D = polyReduc(D,rootD(k));
    end
end