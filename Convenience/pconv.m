function P = pconv(varargin)
% pconv(p1,p2,...) performs vector convolutions for polynomial 
% multiplication
% 
% Graham Holt, May 2026. Updated July 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% pconv(p1,p2,...)
% P = pconv(___)
% 
%% Description
% pconv(p1,p2,...) returns a row vector with the coefficients of the 
% product of input polynomials numeric or symbolic vectors of coefficients 
% using discrete, direct convolution

P = reshape(varargin{1},1,[]);
for k = 2:length(varargin)
    % Steps through all inputs to convolve each vector in input
    n = length(P); m = length(varargin{k});
    M = [P.*reshape(varargin{k},[],1) zeros(m,m-1)];

    P = zeros(1,n+m-1);
    for j = 1:m
        P = P + circshift(M(j,:),j-1);
    end
end