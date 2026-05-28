function [stable, RH] = hurwitz(p)
% hurwitz(p) determines conditions for Hurwitz stability on polynomials
% 
% Graham Holt, March 2026. Updated May 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% hurwitz(p)
% stable = hurwitz(p)
% [stable, RH] = hurwitz(p)
% 
%% Description
% hurwitz(p) takes in a numerical or symbolic vector polynomial "p" and 
% determines if it is Hurwitz stable, or the conditions for such stability.
% If the polynomial "p" is numeric, the output "stable" will be either 0 or
% 1. If not, "stable" will be a vector of symbolic conditions for which "p"
% would be stable
%
% [stable, RH] = hurwitz(p) includes the variable "RH" which contains the
% Routh-Hurwitz table used to determine stability, conditional or otherwise


% Initializes the Routh-Hurwitz Table
n = length(p); m = floor((n+1)/2);
if mod(n,2)==0
    RH = [p(1:2:n); p(2:2:n)];
else
    RH = [p(1:2:n); p(2:2:n) 0];
end
syms epsilon;
RH = sym(RH);
RH(RH(:,1)==0,1) = epsilon;

% Constructs the Routh-Hurwitx Table
while ~all(RH(end,:)==0)
    RH = [RH; zeros(1,m)];
    for j = 1:(m-1)
        RH(end,j) = (1/RH(end-1,1))*det([RH(end-1,1) RH(end-1,j+1);...
                                         RH(end-2,1) RH(end-2,j+1)]);
        RH(end,j) = RH(end,j);
        if RH(end,1)==0 && RH(end-1,2)~=0 && RH(end-2,2)~=0
            RH(end,1) = epsilon;
        end
    end
end
RH(end,:) = [];

% Determines conditions for stability
if ~isa(p,'numeric')
    % Checks the all-positive case
    stable = []; assume(symvar(RH),'real'); assumeAlso(epsilon>0);
    for k = 1:size(RH,1)
        stable = [stable; simplify(RH(k,1)>0)];
        if any(stable==symfalse)
            assume(symvar(RH),'clear');
            break;
        elseif any(stable~=symtrue)
            assumeAlso(stable(stable~=symtrue));
        end
    end
    % Checks the all-negative case, if positive case fails
    if any(stable==symfalse)
        stable = []; assume(symvar(RH),'real'); assumeAlso(epsilon<0);
        for k = 1:size(RH,1)
            stable = [stable; simplify(RH(k,1)<0)];
            if any(stable==symfalse)
                break;
            elseif any(stable~=symtrue)
                assumeAlso(stable(stable~=symtrue));
            end
        end
    end
    if any(stable==symfalse)
        stable = false;
    end
else
    stable = all(sign(RH(:,1))==sign(RH(1,1)));
end

% Removes "epsilon" variable which stands in for zeros
stable = limit(stable,epsilon,0,'right');
RH = limit(RH,epsilon,0,'right');