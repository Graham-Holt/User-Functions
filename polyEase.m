function p = polyEase(x,d,f)
% polyEase(x,d,f) generates a polynomial easing function using derivatives
% at various points
% 
% Graham Holt, April 2026. Updated June 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% polyEase(x,d,f)
% p = polyEase(___)
% 
%% Description
% polyEase(x,d,f) returns the simplest polynomial which has the derivatives 
% "d" equal to "f" at "x"

% Ensures that the polynomial is possible
if any(d >= length(d)-1)
    error(['Order of derivative (',num2str(max(d)),') must be less than order of polynomial (',num2str(length(d)-1),').']);
end
if size(f,1)>1
    f = f.';
end

% Initializes matrix for solution
n = length(d); 
if ~(isa(x,'numeric')&isa(f,'numeric'))
    P = sym(zeros(n));
else
    P = zeros(n);
end

% Constructs matrix and augments with function values
for k = 1:n
    for j = 1:(n-d(k))
        P(k,j+d(k)) = prod(j - 1 + (1:d(k)))*(x(k)^(j-1));
    end
end
Paug = rref([P f.']);

% Checks if the conditions are compatible or redundant
if det(P)==0 && Paug(end,end)==0
    error('One or more conditions are redundant.');
elseif det(P)==0 && Paug(end,end)~=0
    error('Conditions are incompatible.');
end
p = flip(f/(P.'));