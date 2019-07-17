function fnh = dynamicsFunction(expr,states)
% Justin Yim - 2018
% 
% Turns expressions expr in terms of states and time into a function handle
% of time and the state vector;
% INPUTS
% expr: symbolic expressions
% states: cell array of states in order (symbolic variables or characters)
% OUTPUTS
% fnh: function handle taking (t,X) as input

warning('off','symbolic:sym:sym:DeprecateExpressions')
warning('off','symbolic:mupadmex:CallbackError')
warning('off','symbolic:generate:FunctionNotVerifiedToBeValid')

n_states = length(states);
X = sym('X',[1, n_states]);
Xc = cell(1, n_states);
for ii = 1:n_states; Xc{ii} = X(ii); end;

dynState = subs(expr, states, Xc);
dynString1 = char(dynState);
dynString2 = regexprep(dynString1,'X([0-9]+)','X\($1\)');
dynSym = sym(dynString2);
fnh = matlabFunction(dynSym,'Vars',{'t','X'});
