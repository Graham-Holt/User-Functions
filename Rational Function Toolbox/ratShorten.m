function [N,D] = ratShorten(N,D,tol)
% ratShorten(N,D) eliminates leading zeros from rational functions
% 
% Graham Holt, September 2026. Updated September 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% ratShorten(N,D)
% ratShorten(N,D,tol)
% P = ratShorten(___)
% 
%% Description
% ratShorten(N,D) returns the rational function without leading zeros
%
% ratShorten(N,D,tol) allows specification on what to consider a leading
% zero

if ~exist('tol','var')
    tol = 1e-12;
end

N = N(find(abs(N)>=tol,1):end);
D = D(find(abs(D)>=tol,1):end);

if isempty(N) || isempty(D)
    N = 0; D = 0;
elseif isempty(N)
    N = 0; D = 1;
elseif isempty(D)
    N = 1; D = 0;
else
    N = N/D(1); D = D/D(1);
end