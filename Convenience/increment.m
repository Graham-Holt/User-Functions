function v = increment(varargin)
% increment(v,min,max,inc) increments within the provided number system
% 
% Graham Holt, August 2026. Updated August 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% increment(v,min,max,inc)
% increment(v,numSys)
% V = vander2(___)
% 
%% Description
% increment(v,min,max,inc) returns the vector incremented once by "inc"
% between "min" and "max"
%
% increment(v,numSys) allows the min, max, and inc vectors to be stacked.

v = varargin{1};
if nargin == 2
    min = varargin{2}(1,:);
    max = varargin{2}(2,:);
    inc = varargin{2}(3,:);
else
    min = varargin{2};
    max = varargin{3};
    inc = varargin{4};
end

if isscalar(min)
    min = min*ones(1,length(v));
end
if isscalar(max)
    max = max*ones(1,length(v));
end
if isscalar(inc)
    inc = inc*ones(1,length(v));
end

if all(v == max)
    v = min;
    return;
end

v(end) = v(end) + inc(end);

k = 0;
while v(end-k) > max(end-k)
    v(end-k) =  min(end-k);
    v(end-k-1) = v(end-k-1) + inc(end-k-1);
    
    if k+1 == length(v)
        break;
    end
    k = k + 1;
end