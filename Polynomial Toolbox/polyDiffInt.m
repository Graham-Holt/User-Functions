function P = polyDiffInt(p,m,x0,p0)
% polyDiffInt(p,d) performs polynomial differentiation and integration
% 
% Graham Holt, July 2026. Updated July 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% polyDiffInt(p,m)
% polyDiffInt(___,x0,p0)
% P = polyDiffInt(___)
% 
%% Description
% polyDiffInt(p,m) returns a row vector with the coefficients of the m-th
% differ-intergral of the input polynomial
%
% polyDiffInt(___,x0,p0) considers initial/boundary values for the integral
% case with increasing derivative order (p0 = [p(x0) p'(x0) p''(x0) ...])


P = reshape(p,1,[]); n = length(p);
if m>0
    for k = 1:m
        del = (n-k):-1:0;
        P = circshift(P.*del,1);
        P = P(2:end);
    end
else
    if nargin<3
        x0 = zeros(1,-m);
        p0 = zeros(1,-m);
    end
    for k = 1:-m
        int = 1./((n+k-1):-1:1);
        P = P.*int;
        P = [P, p0(1-m-k)-polyval([P,0],x0(1-m-k))];
    end
end