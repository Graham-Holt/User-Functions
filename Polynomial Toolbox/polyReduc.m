function p = polyReduc(p,r)
% polyReduc(p,r) divides out roots from polynomials
% 
% Graham Holt, August 2026. Updated August 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% polyReduc(p,r)
% P = polyReduc(___)
% 
%% Description
% polyReduc(p,r) returns a row vectors with the coefficients of the reduced
% polynomial using synthetic division

m = length(r); P = p;
for k = 1:m
    P = p; n = length(P);
    for j = 2:n
       P(j) = P(j) + P(j-1)*r(k);
    end
    if abs(P(end))>1e-6
        P = p; continue;
    end
    p = P(1:(end-1));
end