function [N,D] = ratComp(varargin)
% ratComp(n1,d1,n2,d2,...) composes rational function in the order 
% provided
% 
% Graham Holt, August 2026. Updated August 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% [N,D] = ratComp(n1,d1,n2,d2,...)
% [N,D] = ratComp(n,d)
% 
%% Description
% ratComp(n1,d1,n2,d2,...) returns a row vector with the coefficients of 
% the composition of the given polynomials in the order r1(r2(...rN(x)...))
%
% ratComp(n,d) allows for "n" and "d" to be a cell array containing 
% polynomial coefficient vectors

if iscell(varargin{1}) && nargin==2
    varargin = collate(varargin{1},varargin{2});
end

N = reshape(varargin{1},1,[]);
D = reshape(varargin{2},1,[]);
for k = 2:(length(varargin)/2)
    n = length(N); m = length(D);
    p = length(varargin{2*k-1});
    q = length(varargin{2*k});
    
    Ntemp = zeros(1,1+(n-1)*(p-1));
    for j = 1:n
        Ntemp = Ntemp + [zeros(1,(j-1)*(p-q)), N(j)*polyMult([repmat(varargin(2*k-1),1,n-j),repmat(varargin(2*k),1,j-1)])];
    end
    N = Ntemp;
    Dtemp = zeros(1,1+(m-1)*(q-1));
    for j = 1:m
        Dtemp = Dtemp + [zeros(1,(j-1)*(p-q)), D(j)*polyMult([repmat(varargin(2*k-1),1,m-j),repmat(varargin(2*k),1,j-1)])];
    end
    D = Dtemp;

    if n>m
        D = polyMult([{D},repmat(varargin{2*k},1,n-m)]);
    elseif n<m
        N = polyMult([{N},repmat(varargin{2*k},1,m-n)]);
    end
end