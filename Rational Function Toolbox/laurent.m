function [L,radiusConv] = laurent(N,D,nTerms,z0)
% laurent(N,D) computes principle Laurent series about its poles
% 
% Graham Holt, September 2026. Updated September 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% laurent(N,D)
% laurent(___,nTerms)
% laurent(___,z0)
% L = laurent(___)
% [L,radiusConv] = laurent(___)
% 
%% Description
% laurent(N,D) returns a cell column with 6 terms of the principle Laurent
% series about each of its poles (including the 0th term)
%
% laurent(___,nTerms) returns the specified number of terms, or fewer
%
% laurent(___,z0) returns a cell column of the principle Laurent series
% about each value in "z0"
%
% [L,radiusConv] = laurent(___) returns the radius of convergence for each
% of the series

if ~exist('nTerms','var')
    nTerms = 6;
end
if ~exist('z0','var')
    z0 = unique(polyRoots(D));
end

for k = 1:length(z0)    
    if z0(k) == 0
        W = {1, [1 0]};
    else
        W = {[z0(k) 1], [1 0]};
    end
    
    [Nw,Dw] = ratComp(N,D,W{1},W{2});
    
    L{k,:}(1) = ratval(Nw,Dw,0);
    for j = 2:nTerms
        [Nw,Dw] = ratDiff(Nw,Dw*(j-1),1);
    
        L{k}(j) = ratval(Nw,Dw,0);
    end
    L{k} = L{k}(~isnan(L{k}));
    
    denomRoots = polyRoots(D); potentialRadii = abs(denomRoots - z0(k));
    potentialRadii = potentialRadii(potentialRadii > 1e-12);
    radiusConv(k,:) = min(potentialRadii); 
    if isempty(radiusConv(k)) 
        radiusConv(k) = Inf;
    end
end