function g=gradientfd(a,test_function)
% Objective: Generates Gradient of onjective function at specific point
%-----------------------------------------------------------------------
% f=LagGrad(a,test_function)
% where a=input vector
%       test_function=objective function
%-----------------------------------------------------------------------
% Output: f= gradient column vector
%-----------------------------------------------------------------------

% Code by:
% Salil Sharma
% For the project implementation in IE 538 course
% Spring 2017
%-----------------------------------------------------------------------

ep=0.0001; % step size for numerical differentiation
l=length(a); % length of input vector
val=test_function(a); % value of obj function
for i=1:l
    x1=a;
    x1(i)=a(i)+ep;
    v=test_function(x1);
    df(i)=(v-val)/ep; %forward difference
end 
g=df; % returns gradient vector
end
