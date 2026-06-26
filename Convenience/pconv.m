function P = pconv(p1,p2)
% pconv(p1,p2) performs vector convolutions
% 
% Graham Holt, May 2026. Updated June 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% pconv(p1,p2)
% P = pconv(p1,p2)
% 
%% Description
% pconv(p1,p2) takes two numeric or symbolic vectors of polynomial 
% coefficients and returns a row vector with the coefficients of their 
% product by discrete, direct convolution

% Checks and re-orients inputs to match algorithm
if ~(isvector(p1) && isvector(p2))
    error('Inputs must be vectors.');
end

p1 = reshape(p1,[],1); n1 = length(p1);
p2 = reshape(p2,1,[]); n2 = length(p2);

% Builds table of products and combines like-terms
M = [p1.*p2 zeros(n1,n1-1)];
P = zeros(1,n1+n2-1);
for k = 1:n1
    P = P + circshift(M(k,:),k-1);
end