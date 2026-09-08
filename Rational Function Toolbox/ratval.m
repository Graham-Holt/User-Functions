function R = ratval(N,D,x)
% ratval(N,D,x) evaluates ratiopnal functions
% 
% Graham Holt, September 2026. Updated September 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% ratval(N,D,x)
% P = ratval(___)
% 
%% Description
% ratval(N,D,x) returns a vector of the rational function N/D evaluated at 
% the points in "x"

R = polyval(N,x)./polyval(D,x);