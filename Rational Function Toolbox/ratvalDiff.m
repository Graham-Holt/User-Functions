function R = ratvalDiff(N,D,x,d)
% ratvalDiff(N,D,x) evaluates rational functions and their derivatives.
% 
% Graham Holt, November 2026. Updated November 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% ratvalDiff(N,D,x)
% ratvalDiff(___,d)
% P = ratvalDiff(___)
% 
%% Description
% ratvalDiff(N,D,x) returns the same result as the usual ratval() function
%
% ratvalDiff(N,D,x,d) returns a similar result with the value of the (d)th
% derivative of the rational function at (x)

x = reshape(x,[],1);

if nargin < 3
    d = 0;
end

for k = 1:length(d)
    R(:,k) = ratval(ratDiff(N,D,d(k)),x);
end