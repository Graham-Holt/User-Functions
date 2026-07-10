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
% polyInterp(x,d,f) returns the simplest polynomial which has the 
% derivatives "d" equal to "f" at "x"

% Ensures that there exists an interpolating solution (no least-squares)
if ~(length(d)==length(x) && length(f)==length(x))
    error('Vetcor inputs must have the same length');
end
if any(d>=length(d)) || any(d~=round(d))
    error('Order of derivative must be an integer less than or equal to the order of polynomial.');
end
x = reshape(x,[],1); f = reshape(f,[],1);

% Initializes matrix for solution
n = length(x); P = zeros(n);
if ~(isa(x,'numeric')&isa(f,'numeric'))
    P = sym(P);
end

% Constructs matrix and augments with function values
del = n - (1:n);
for k = 1:n
    % Performs deriviative on row of Vandemonde matrix
    P(k,:) = x(k).^del;
    for j = 1:d(k)
       P(k,:) = circshift(del.*P(k,:),1);
    end
    
    % Corrects for vanishing effect when x == 0
    if x(k)==0
        P(k,:) = [zeros(1,n-1) 1];
    end

    % Shifts Vandemonde matrix to align with coefficients
    P(k,:) = circshift(P(k,:),-d(k));
end
Paug = rref([P f]); r = rank(P);

% Finds minimum norm solution
if r==n 
    p = Paug(:,end);
elseif r~=n && Paug(end,end)==0
    P = Paug(1:r,1:(end-1));
    f = Paug(1:r,end);
    p = P.'/(P*P.')*f;
else
    error('One or more conditions are incompatible.');
end

% Reduces polynomial to minimal degree
p = p(find(p~=0,1):end).';