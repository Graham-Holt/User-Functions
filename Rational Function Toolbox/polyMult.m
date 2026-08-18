function P = polyMult(varargin)
% polyMult(p1,p2,...) performs vector convolutions for polynomial 
% multiplication
% 
% Graham Holt, May 2026. Updated August 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% polyMult(p1,p2,...)
% polyMult(p)
% P = polyMult(___)
% 
%% Description
% polyMult(p1,p2,...) returns a row vector with the coefficients of the 
% product of input polynomials numeric or symbolic vectors of coefficients 
% using discrete, direct convolution
%
% polyMult(p) allows for "p" to be a cell array containing polynomial
% coefficient vectors

if iscell(varargin{1}) && nargin==1
    varargin = reshape(varargin{1},1,[]);
end

% Steps through all inputs to convolve each vector in input
P = reshape(varargin{1},1,[]);
for k = 2:length(varargin)
    n = length(P); m = length(varargin{k});

    M = [P.*reshape(varargin{k},[],1) zeros(m,m-1)];

    P = zeros(1,n+m-1);
    for j = 1:m
        P = P + circshift(M(j,:),j-1);
    end
end