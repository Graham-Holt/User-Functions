function [funcUnc, funcEvalX] = propUnc(func,xName,xValue,xUnc,constName,constValue,method)
% propUnc(func,xName,xValue,xUnc) evaluates the uncertainty and value of a 
% function near a set of inputs with given uncertainties
% 
% Graham Holt, March 2025. Updated July 2026
% Embry-Riddle Aeronautical University
% 
%% Syntax
% propUnc(func,xName,xValue,xUnc)
% propUnc(___,constName,constValue)
% propUnc(___,method)
% funcUncertainty = propUnc(___)
% [funcUncertainty, funcEvalX] = propUnc(___)
% 
%% Description
% ****All inputs are either numeric arrays or a row vector of strings.****
%
% propUnc(func,xName,xValue,xUnc) takes the function "func" with variables 
% "xName", and returns the uncertainty of the function with variable values 
% "xValue" and uncertainties "xUnc".
%
% propUnc(___,constName,constValue) clears up some notation by defining 
% constants used in the function "func" which don't need different values. 
% If not included, use empty string or cell array.
%
% propUnc(___,method) sets whether or not the output uses the standard,
% symmetric RMS method using partial derivatives (default), or computes the
% asymmetric upper and lower bounds of the values ('abs').
%
% [funcUnc, funcEvalX] = propUnc(___) includes the output of the function 
% values for the given inputs.
%
%% Example
% Consider a calculation of air density from pressure, temperature, and the
% specific gas constant.
% 
%  propUnc("p/(R*T)",["p","T","R"],...
%          [pres, temp, 287*ones(length(pres),1)],...
%          [pres_Unc,temp_Unc,zeros(length(pres),1])
%
% We can tidy this up by passing "R" as a constant value.
%
%  propUnc("p/(R*T)",["p","T"],[pres, temp],[pres_Unc, temp_Unc],"R",287)
%
% If we had upper and lower bounds for pressure and temperature, we can use
% the following.
%
%  propUnc("p/(R*T)",["p","T"],[pres, temp],...
%          [pres_Unc_lo, pres_Unc_hi, temp_Unc_lo, temp_Unc_hi],...
%          "R",287)
%
% Additionally, by using the 'abs' method, we can get asymmetric upper and 
% lower bounds for the resulting air density.

% Converts strings into usable symbolic objects
xName = str2sym(xName); func = str2sym(func);

% Pre-inserts constants where possible
if exist('constName','var') && exist('constValue','var')
    constName = str2sym(constName);
    func = subs(func,constName,constValue);
end

% Selects default method
symmetric = size(xValue,2)==size(xUnc,2);
if ~exist('method','var')
    method = 'rms';
end

% Evaluates function at the given inputs
for k = 1:size(xValue,1)
    funcEvalX(k,:) = double(subs(func,xName,xValue(k,:)));
end

if strcmpi(method,'abs')
% Evaluates at all possible uncertainies for comparision

% Puts symmetric uncertainies into asymmetric form
if symmetric
    collate = reshape([0; size(xValue,2)]+(1:size(xValue,2)),1,[]);
    xUnc = [-xUnc xUnc]; xUnc = xUnc(:,collate);
end

for k = 1:(2^size(xValue,2))    
    % Selects each combination of uncertainties to comparison
    selectXUncertainty = (dec2bin(k-1,size(xValue,2))-'0') - 1 + (2:2:size(xUnc,2));
    for j = 1:size(xValue,1)
        funcEvalX_perturb(j,k) = double(subs(func,xName,xValue(j,:)+xUnc(j,selectXUncertainty)));
    end
end
funcUnc = double([min(funcEvalX_perturb,[],2) max(funcEvalX_perturb,[],2)] - funcEvalX);

if symmetric
    funcUnc = max(abs(funcUnc),[],2);
end

else
% Uses partial derivatives to approximate effect of uncertainities

% Considers largest deviation for uncertainty
if ~symmetric
    xUnc = cat(3,xUnc(:,1:2:end),xUnc(:,2:2:end));
end

for k = 1:length(xName)
    for j = 1:size(xValue,1)
        funcUnc(j,k,:) = double(subs(diff(func,1,xName(k)),xName,xValue(j,:))).*xUnc(j,k,:);
    end
end
funcUnc = sqrt(sum(funcUnc.^2,2));
funcUnc = reshape(funcUnc.*sign(xUnc(1,1,:)),size(funcUnc,1),[],1);

end