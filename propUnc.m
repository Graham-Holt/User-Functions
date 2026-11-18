function [funcUnc, funcEvalX] = propUnc(func,xValue,xUnc,symmetric)
% propUnc(func,xValue,xUnc) evaluates the uncertainty and value of a 
% function near a set of inputs with given uncertainties
% 
% Graham Holt, March 2025. Updated November 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% propUnc(func,xValue,xUnc)
% propUnc(___,symmetric)
% funcUncertainty = propUnc(___)
% [funcUnc, funcEvalX] = propUnc(___)
% 
%% Description
% propUnc(func,xValue,xUnc) takes the function "func" with variables, and 
% returns the uncertainty of the function with variable values "xValue" and 
% uncertainties "xUnc".
%
% propUnc(___,symmetric) sets whether or not the output returns a singular,
% maximal value for uncertainty or asymmetric upper and lower bounds.
%
% [funcUnc, funcEvalX] = propUnc(___) includes the output of the function 
% values for the given inputs.
%
%% Example
% Consider a calculation of air density from pressure, temperature, and the
% specific gas constant.
% 
%  propUnc(@(x) x(1)./(R*x(2)), [pres, temp],[pres_Unc,temp_Unc])
%
% If we had upper and lower bounds for pressure and temperature, we can use
% the following.
%
%  propUnc(@(x) x(1)./(R*x(2)),[pres, temp],...
%          [pres_Unc_lo, pres_Unc_hi, temp_Unc_lo, temp_Unc_hi])

% Selects defaults and pre-inserts constants where possible
if nargin < 4
    symmetric = size(xValue,2)==size(xUnc,2);
end

% Puts symmetric uncertainies into asymmetric form
if size(xValue,2) == size(xUnc,2)
    xUnc = collate(-xUnc,xUnc);
end
if isvector(xUnc)
    xUnc = repmat(xUnc,[size(xValue,1) 1]);
end

% Evaluates at all possible uncertainies for comparision
funcEvalX = func(xValue);

for k = 1:(2^size(xValue,2))    
    % Selects each combination of uncertainties to comparison
    selectXUncertainty = (dec2bin(k-1,size(xValue,2))-'0') - 1 + (2:2:size(xUnc,2));
    for j = 1:size(xValue,1)
        funcEvalX_perturb(j,k) = func(xValue(j,:) + xUnc(j,selectXUncertainty));
    end
end
funcUnc = [min(funcEvalX_perturb,[],2) max(funcEvalX_perturb,[],2)] - funcEvalX;

if symmetric
    funcUnc = max(abs(funcUnc),[],2);
end

end

function V = collate(varargin)
% collate(v1,v2,...) collates same-size vectors.
% 
% Graham Holt, August 2026. Updated August 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% collate(v1,v2,...)
% collate(A1,A2,...)
% V = collate(___)
% 
%% Description
% collate(v1,v2,...) returns a row vector of the collated elements of the
% input vectors
%
% collate(A1,A2,...) takes in arrays as input and returns a new array with
% the columns collated

n = length(varargin); m = size(varargin{1},2);
for k = 1:m
    for j = 1:n
        V(:,n*(k-1) + j) = varargin{j}(:,k);
    end
end

end