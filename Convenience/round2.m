function y = round2(varargin)
% round2(x,roundTo) rounds to the nearest decimal
% 
% Graham Holt, August 2026. Updated August 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% round2(x,roundTo)
% round2(___,type)
% round2(___,dir)
% y = round2(___)
% 
%% Description
% round2(x,roundTo) rounds the input to the specified number of digits
%
% round2(___,type) distinguishes what's being rounded to:
%   0 - round to number of digits after the decimal (default)
%   1 - round to number of significant figures
%   2 - round to adjacent multiple of "roundTo"
%
% round2(___,dir) specifies which direction to round as a number 
% (Inf -> ceil, -Inf -> floor, 0 -> fix, etc.)

if ~isnumeric([varargin{:}])
    error('Inputs must be numeric.');
end
if nargin<2
    error('Too few input arguments.');
end
if nargin>4
    error('Too many input arguments.');
end

x = varargin{1}; 
roundTo = varargin{2};

if nargin<3
    type = 0;
else
    type = varargin{3};
end

if nargin<4
    signDir = 1 - 2*(2*mod(x,roundTo) < roundTo);
    dir = signDir*Inf;
else
    dir = varargin{4};
end

if type~=2 || isempty(type)
    roundTo = 10^(type*(floor(log10(x)) + 1) - roundTo);
end

signRound = 1 - 2*(x < dir);
y = x - signRound*mod(signRound*x,roundTo);