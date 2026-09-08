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

if ~exist('tol','var')
    tol = 1e-6;
end

rootD = polyRoots(D,tol);

for k = 1:length(rootD)
    n = length(N);

    N = polyReduc(N,rootD(k));
    if length(N)~=n
        D = polyReduc(D,rootD(k));
    end
end

N = N(find(N~=0,1):end);
D = D(find(D~=0,1):end);
if isempty(N) || isempty(D)
    N = 0; D = 0;
elseif isempty(N)
    N = 0; D = 1;
elseif isempty(D)
    N = 1; D = 0;
end