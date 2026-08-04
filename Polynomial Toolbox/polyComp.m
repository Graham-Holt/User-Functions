function P = polyComp(varargin)
% polyComp(p1,p2,...) composes polynomials in the order provided.
% 
% Graham Holt, July 2026. Updated July 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% polyComp(p1,p2,...)
% polyComp(p)
% P = polyComp(___)
% 
%% Description
% polyComp(p1,p2,...) returns a row vector with the coefficients of the 
% composition of the given polynomials in the order p1(p2(...pN(x)...))
%
% polyComp(p) allows for "p" to be a cell array containing polynomial
% coefficient vectors

if iscell(varargin{1}) && nargin==1
    varargin = reshape(varargin{1},1,[]);
end

P = reshape(varargin{1},1,[]);
for k = 2:length(varargin)
    n = length(P); m = length(varargin{k});
    M = zeros(1,1 + (n-1)*(m-1));
    for j = 1:(n-1)
        M = M + [zeros(1, (j-1)*(m-1)), P(j)*polyConv(repmat(varargin(k),1,n-j))];
    end
    P = M + [zeros(1, (n-1)*(m-1)) P(n)];
end