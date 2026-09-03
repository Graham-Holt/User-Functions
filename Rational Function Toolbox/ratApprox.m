function [N,D] = ratApprox(func,interval,order)
% ratApprox(x,f) generates rational approximations using Chebeshev nodes
% 
% Graham Holt, August 2026. Updated September 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% ratApprox(func,interval,order)
% [N,D] = ratApprox(___)
% 
%% Description
% ratApprox(func,interval,order) returns the simplest polynomial which
% interpolate the function values within the given interval to the desired
% order.
%
% [N,D] = ratApprox(___) returns the simplest rational function when the
% order is a vector [numeratorOrder, denominatorOrder]

midInt = (interval(2) + interval(1))/2;
radiusInt = abs(interval(2) - interval(1))/2;

order = [reshape(order,1,[]) zeros(1,2-length(order))];
numPoints = order(1) + order(2) + 1;

if numPoints==1
    x_est = midInt;
else
    x_est = midInt + radiusInt*cos((0:(numPoints-1))*pi/(numPoints-1)).';
end

[N,D] = ratInterp(x_est,func(x_est),order(1)-order(2));