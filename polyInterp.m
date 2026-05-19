function p = polyInterp(x,d,f)
% polyInterp(x,d,f) generates interpolating polynomials using derivatives
% at various points
% 
% Graham Holt, April 2026. Updated May 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% polyInterp(x,d,f)
% p = polyInterp(___)
% 
%% Description
% polyInterp(x,d,f) takes returns the simplest polynomial which has the 
% derivatives "d" equal to "f" at "x"

if any(d >= length(d)-1)
    error(['Order of derivative (',num2str(max(d)),') must be less than order of polynomial (',num2str(length(d)-1),').']);
end

n = length(d); P = zeros(n);
for k = 1:n
    for j = 1:(n-d(k))
        P(k,j+d(k)) = prod(j - 1 + (1:d(k)))*(x(k)^(j-1));
    end
end
Paug = rref([P f']);
if det(P)==0 && Paug(end,end)==0
    error('One or more conditions are redundant.');
elseif det(P)==0 && Paug(end,end)~=0
    error('Conditions are incompatible.');
end
p = flip(f/(P'));