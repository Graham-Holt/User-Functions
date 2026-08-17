function y = roundd(x,roundTo,type,dir)
% roundd(x,roundTo) rounds to the nearest decimal
% 
% Graham Holt, August 2026. Updated August 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% roundd(x,roundTo)
% roundd(___,type)
% roundd(___,dir)
% y = roundd(___)
% 
%% Description
% roundd(x,roundTo) rounds the input to the specified number of digits
%
% roundd(___,type) distinguishes rounding to "roundTo" digits (0) or the
% nearest multiple of "roundTo" (1)
%
% roundd(___,dir) specifies which direction to round as a number 
% (Inf -> ceil, -Inf -> floor, 0 -> fix, etc.)

if ~exist('type','var')
    type = 0;
end
if ~type
    roundTo = 10^(-roundTo);
end

if nargin<4
    if 2*mod(x,roundTo) > roundTo
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

y = x - s*mod(s*x,roundTo);