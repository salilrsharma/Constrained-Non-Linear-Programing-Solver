function g=LagGrad(a,confunction)
% Objective: Generates Gradient of constraint functions at specific point
%-----------------------------------------------------------------------
% f=LagGrad(a,confunction)
% where a=input vector
%       confunction=constraints arrnaged in a nx1 matrix
%-----------------------------------------------------------------------
% Output: f= Column array of gradients for individual constraint 
%            functions arranged
%-----------------------------------------------------------------------

% Code by:
% Salil Sharma
% For the project implementation in IE 538 course
% Spring 2017
%-----------------------------------------------------------------------

ep=0.0001; % Step size for numerical differentiation
for i=1:length(a)
    x1=a;
    x1(i)=a(i)+ep;
    b(i,:)=(confunction(x1)-confunction(a))/ep; % Forward difference method
end
g=b; % return gradients in individual columns of one matrix 
end
