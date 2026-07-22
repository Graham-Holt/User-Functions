function [stable, tabRH] = hurwitz(p)
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
% [stable, tabRH] = hurwitz(p) includes the variable "tabRH" which contains 
% the Routh-Hurwitz table used to determine stability, conditional or 
% otherwise


% Keeps previous assumptions in tact
prevAssume = assumptions;

% Removes need to check both positive and negative case
p = sign(p(1))*p;

% Initializes the Routh-Hurwitz Table
polyLen = length(p); tabWidth = floor((polyLen+1)/2);
tabRH = sym(zeros(polyLen,tabWidth));
if mod(polyLen,2)==0
    tabRH(1:2,:) = [p(1:2:polyLen); p(2:2:polyLen)];
else
    tabRH(1:2,:) = [p(1:2:polyLen); p(2:2:polyLen) 0];
end
syms epsilon; assume(epsilon>0);
if tabRH(2,1)==0  tabRH(2,1) = epsilon; end

% Constructs the Routh-Hurwitx Table
for k = 3:polyLen
    for j = 1:(tabWidth-1)
        tabRH(k,j) = (1/tabRH(k-1,1))*det([tabRH(k-1,1) tabRH(k-1,j+1);...
                                     tabRH(k-2,1) tabRH(k-2,j+1)]);
        if k==polyLen && tabRH(k,1)==0
            stable = false;
            warning('Marginal stability detected.');
            return;
        elseif tabRH(k,1)==0
            tabRH(k,1) = epsilon;
        end
    end
end

if ~isa(p,'numeric')
    
    % Initializes conditions and Routh-Hurwitz table
    tabRH = simplify(tabRH);
    stable = sym(zeros(size(tabRH,1),1));
    
    % Determines symbolic conditions for stability
    for k = 1:size(tabRH,1)
        stable(k) = simplify(tabRH(k,1)>0);
        if stable(k)==symfalse
            stable = false;
            break;
        end
    end

else
    % Determines stability for numeric case
    stable = all(sign(tabRH(:,1))==sign(tabRH(1,1)));
end

% Removes "epsilon" variable, which stands in for zeros
stable = limit(stable,epsilon,0,'right');
tabRH = limit(tabRH,epsilon,0,'right');
assume(epsilon,'clear');

if ~isa(p,'numeric')

% Simplifies symbolic conditions for stability
for k = 1:length(stable)
    if stable(k)==symtrue
        continue;
    end
    assume([prevAssume stable(k)]);
    stable([1:(k-1) (k+1):end]) = simplify(stable([1:(k-1) (k+1):end]));
end

end

% Reinstates previously defined assumptions
assume(prevAssume);