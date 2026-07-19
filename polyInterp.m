function p = polyInterp(x,d,f)
% polyInterp(x,d,f) generates the simplest polynomial function defined by 
% its derivatives at various points
% 
% Graham Holt, April 2026. Updated July 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% polyInterp(x,d,f)
% p = polyInterp(___)
% 
%% Description
% polyInterp(x,d,f) returns the simplest polynomial which has derivatives 
% "d" equal to "f" at "x"

% Ensures that there exists an interpolating solution (no least-squares)
if ~(length(d)==length(x) && length(f)==length(x))
    error('Vetcor inputs must have the same length');
end
if any(d>=length(d)) || any(d~=round(d))
    error('Order of derivative must be an integer less than or equal to the order of polynomial.');
end
x = reshape(x,[],1); f = reshape(f,[],1);

% Initializes Vandermonde matrix
n = length(x); P = zeros(n);
if ~(isa(x,'numeric')&isa(f,'numeric'))
    P = sym(P);
end

% Constructs matrix and augments with function values
del = n - (1:n);
for k = 1:n
    % Create shifted Vandemonde matrix
    P(k,:) = [x(k).^del((1+d(k)):n) zeros(1,d(k))];

    % Multiply coefficients using power rule
    for j = 1:d(k)
        P(k,:) = P(k,:).*circshift(del,1-j);
    end
end
Paug = rref([P f]); r = rank(P);

if any(Paug((r+1):n,n+1)~=0)
    error(['(',num2str(sum(Paug((r+1):n,n+1)~=0)),') conditions are incompatible.']);
end

% Finds minimum norm solution
P = Paug(1:r,1:n);
f = Paug(1:r,n+1);
p = P.'/(P*P.')*f;

% Reduces polynomial to minimal degree
p = p(find(p~=0,1):end).';