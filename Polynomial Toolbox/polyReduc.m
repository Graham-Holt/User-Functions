function p = polyReduc(p,r,tol)
% polyReduc(p,r) divides out roots from polynomials
% 
% Graham Holt, August 2026. Updated August 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% polyReduc(p,r)
% polyReduc(___,tol)
% P = polyReduc(___)
% 
%% Description
% polyReduc(p,r) returns a row vectors with the coefficients of the reduced
% polynomial using synthetic division
%
% polyReduc(___,tol) tunes the tolerance for an input root to be considered
% a root of the given polynomial

if ~exists('tol','var')
    tol = 1e-6;
end

m = length(r); P = p;
for k = 1:m
    P = p; n = length(P);
    if abs(polyval(P,r(k)))>tol
        continue;
    end
    for j = 2:n
       P(j) = P(j) + P(j-1)*r(k);
    end
    p = P(1:(end-1));
end