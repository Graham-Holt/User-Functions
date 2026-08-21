function [r,isFullRank] = rankFF(A,modulo)
% rankFF(A,modulo) computes matrix rank in finite fields
% 
% Graham Holt, August 2026. Updated August 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% rankFF(A,modulo)
% r = rankFF(___)
% [r,isFullRank] = rankFF(___)
% 
%% Description
% rankFF(A,modulo) returns the rank of the matrix "A" in the finite field
% "modulo"
%
% [r,isFullRank] = rankFF(___) additionally returns a boolean that is true
% when the matrix "A" is full rank in the finite field "modulo"

r = size(A,1) - sum(ismember(rrefFF(A,modulo),zeros(1,size(A,2)),'rows'));

isFullRank = r == min(size(A));