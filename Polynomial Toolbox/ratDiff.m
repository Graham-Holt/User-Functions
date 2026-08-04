function [N,D] = ratDiff(n,d,m)
% ratDiff(n,d,m) performs differentiation on rational functions
% 
% Graham Holt, August 2026. Updated August 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% ratDiff(n,d,m)
% P = ratDiff(___)
% 
%% Description
% ratDiff(n,d,m) returns two row vectors with the coefficients of the 
% m-th derivative of the input rational

if m<=0
    error('Derivative order must be a positive integer');
end

N = n; D = d;
for k = 1:m
    N = polyConv(polyDiffInt(N,1),D) - polyConv(N,polyDiffInt(D,1));
    D = polyConv(D,D);
end

