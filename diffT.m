function diffed = diffT(expr,vars,num)
% Justin Yim - 2018
% 
% Differentiate expression expr in terms of time-dependent variables vars
% with respect to time num times.
% 
% INPUTS:
% expr: symbolic expression
% vars: cell array of chars or symbolic variables
% num: integer number of times to differentiate
% OUTPUTS:
% diffed: symbolic expression with time derivatives written as
%   [char(vars{i}),'d'] for i = 1:length(vars)

warning('off','symbolic:sym:sym:DeprecateExpressions')

syms t

% Write variables as functions of time
for ii = 1:length(vars)
    expr = subs(expr,vars{ii},sym([char(vars{ii}),'(t)']));
end

% Differentiate
dexpr = expr;
for ii = 1:num
    dexpr = diff(dexpr,t);
end

% Replace with correct variable names
dexpr_str = char(dexpr);
for ii = num:-1:1
    find_str = ['diff\((\w+)\(t\)',repmat(', t',[1,ii]),'\)'];
    repl_str = ['$1',repmat('d',[1,ii])];
    dexpr_str = regexprep(dexpr_str,find_str,repl_str);
end
dexpr_str = regexprep(dexpr_str,'\(t\)',''); % Remove explicit time dependence

diffed = sym(dexpr_str);