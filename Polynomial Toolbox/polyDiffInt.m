function P = polyDiffInt(p,d,x0,p0)
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

% Ensures that inputs are valid and derivative is non-trivial
if nargin<3 && d<0
    x0 = zeros(1,-d);
    p0 = zeros(1,-d);
end
if ~(length(x0)==length(p0))
    error('x0 and p0 must be the same length.');
end
if d~=round(d)
    error('Order of differ-integral must be an integer.');
end
if d>=length(p)
    P = 0;
    return;
end

P = reshape(p,1,[]); n = length(p);
if d>0
    for k = 1:d
        del = (n-k):-1:1;
        P = P(1:(end-1)).*del;
    end
else
    for k = 1:-d
        int = 1./((n+k-1):-1:1);
        P = [P.*int, p0(1-d-k)-polyval([P.*int,0],x0(1-d-k))];
    end
end