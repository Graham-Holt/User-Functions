function [stable, RH] = hurwitz(p)
% hurwitz(p) determines conditions for Hurwitz stability on polynomials
% 
% Graham Holt, March 2026. Updated June 2026
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


% Keeps previous assumptions in tact
prevAssume = assumptions;

% Initializes the Routh-Hurwitz Table
n = length(p); m = floor((n+1)/2);
if mod(n,2)==0
    RH = [p(1:2:n); p(2:2:n)];
else
    RH = [p(1:2:n); p(2:2:n) 0];
end
syms epsilon;
RH = sym(RH); RH(RH(:,1)==0,1) = epsilon;

% Constructs the Routh-Hurwitx Table
for k = 3:length(p)
    RH = [RH; zeros(1,m)];
    for j = 1:(m-1)
        RH(k,j) = (1/RH(k-1,1))*det([RH(k-1,1) RH(k-1,j+1);...
                                     RH(k-2,1) RH(k-2,j+1)]);
        if k==length(p) && all(RH(k,:)==0)
            stable = false;
            warning('Marginal stability detected.');
            return;
        elseif RH(k,1)==0
            RH(k,1) = epsilon;
        end
    end
end
RH(end,:) = [];

% Determines conditions for stability
if ~isa(p,'numeric')
    RH = simplify(RH);
    % Checks the all-positive case
    stable = sym(zeros(size(RH,1),1)); assume(epsilon>0);
    for k = 1:size(RH,1)
        stable(k) = simplify(RH(k,1)>0);
        if stable(k)==symfalse
            assume(epsilon,'clear');
            break;
        end
    end
    % Checks the all-negative case, if positive case fails
    if any(stable==symfalse)
    stable = sym(zeros(size(RH,1),1)); assume(epsilon<0);
    for k = 1:size(RH,1)
        stable(k) = simplify(RH(k,1)<0);
        if stable(k)==symfalse
            assume(epsilon,'clear');
            break;
        end
    end
    end
    if stable(k)==symfalse
        stable = false;
    end
else
    assume(epsilon>0);
    stable = all(sign(RH(:,1))==sign(RH(1,1)));
end

% Removes "epsilon" variable which stands in for zeros
stable = limit(stable,epsilon,0,'right');
RH = limit(RH,epsilon,0,'right');
assume(epsilon,'clear');

if ~isa(p,'numeric')
for k = 1:length(stable)
    if stable(k)==symtrue
        continue;
    end
    assume([prevAssume stable(k)]);
    stable([1:(k-1) (k+1):end]) = simplify(stable([1:(k-1) (k+1):end]));
end
end
assume(prevAssume);