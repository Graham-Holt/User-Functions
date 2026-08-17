function y = roundTo(x,rounding,type,dir)
% roundTo(x,rounding) rounds to the nearest decimal
% 
% Graham Holt, August 2026. Updated August 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% roundTo(x,rounding)
% roundTo(___,type)
% roundTo(___,dir)
% y = roundTo(___)
% 
%% Description
% roundTo(x,rounding) rounds the input to the specified number of digits
%
% roundTo(___,type) distinguishes what's being rounded to:
%   0 - round to number of digits after the decimal (default)
%   1 - round to number of significant figures
%   2 - round to the nearest multiple of "rounding"
%
% roundTo(___,dir) specifies which direction to round as a number 
% (Inf -> ceil, -Inf -> floor, 0 -> fix, etc.)

if ~exist('type','var')
    type = 0;
end
if type~=2
    rounding = 10^(type*floor(log10(x)) - rounding);
end

if nargin<4
    if 2*mod(x,rounding) >= rounding
        dir = Inf;
    else
        dir = -Inf;
    end
end
    
if x > dir
    s = 1;
else
    s = -1;
end

y = x - s*mod(s*x,rounding);