function P = pconv(p1,p2)
% pconv(p1,p2) performs vector convolutions
% 
% Graham Holt, May 2026. Updated May 2026
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

if ~(isvector(p1)&&isvector(p2))
    error('Inputs must be vectors.');
end
if xor(size(p1,1)>1,size(p2,1)>1)
    p1 = p1.';
end

P = p1.*p2.';
P = [P zeros(size(P,1),size(P,1)-1)];

for k = 1:size(P,1)
    P(k,:) = circshift(P(k,:),k-1);
end
P = sum(P,1);