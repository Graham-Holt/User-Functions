function P = polyShorten(p,tol)
% polyShorten(p) eliminates leading zeros from polynomials
% 
% Graham Holt, September 2026. Updated September 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% polyShorten(p)
% polyShorten(p,tol)
% P = polyShorten(___)
% 
%% Description
% polyShorten(p) returns the polynomial without leading zeros
%
% polyShorten(p,tol) allows specification on what to consider a leading
% zero

if ~exist('tol','var')
    tol = 1e-12;
end

P = p(find(abs(p)>=tol,1):end);
if isempty(P)
    P = 0;
end