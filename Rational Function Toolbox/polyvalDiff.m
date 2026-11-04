function P = polyvalDiff(p,x,d)
% polvalDiff(p,x) evaluates polynomials and their derivatives.
% 
% Graham Holt, October 2026. Updated October 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% polvalDiff(p,x)
% polvalDiff(p,x,d)
% P = polvalDiff(___)
% 
%% Description
% polvalDiff(p,x) returns the same result as the usual polyval() function
%
% polvalDiff(p,x,d) returns a similar result with the value of the (d)th
% derivative of the polynomial at (x)

if nargin < 3
    d = zeros(1,length(x));
end
if length(d) ~= length(x)
    error('(d) must have the same length as (x)');
end

for k = 1:length(x)
    P(k,:) = polyval(polyDiffInt(p,d(k)),x(k));
end