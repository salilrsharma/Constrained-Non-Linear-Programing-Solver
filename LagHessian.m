function f=LagHessian(a,confunction)
% Objective: Generates Hessian of constraint functions at specific point
%-----------------------------------------------------------------------
% f=LagHessian(a,confunction)
% where a=input vector
%       confunction=constraints arrnaged in a nx1 matrix
%-----------------------------------------------------------------------
% Output: f= Struct array of Hesian matrices for individual constraint 
%            functions arranged
%-----------------------------------------------------------------------

% Code by:
% Salil Sharma
% For the project implementation in IE 538 course
% Spring 2017
%-----------------------------------------------------------------------

l=length(a); %Hessian would be lxl matrix
ep=0.0001; % Step size for numerical differentiation
valf=confunction(a); % value of constraint function at a point 'a'
k=length(valf); %k denotes no of Hessian to be produced on constraints
ep2=ep*ep; 
ep3=4*ep*ep;
for i=1:length(a)
    x1=a;
    x1(i)=a(i)-ep;
    x2=a;
    x2(i)=a(i)+ep;
    z1=(confunction(x2)-2*valf+confunction(x1))/ep2;
    for m=1:k
        s(m).hf(i,i)=z1(m); % Struct array to store the diagonal values of 
                            % hessians of individual contraint functions
    end
    j=i+1;
    while j<=length(a) % to compute rest of the elements in hessian
        x1(j)=a(j)-ep;
        x2(j)=a(j)+ep;
        v4=confunction(x1);
        v1=confunction(x2);
        x1(j)=x1(j)+2*ep;
        x2(j)=x2(j)-2*ep;
        v2=confunction(x1);
        v3=confunction(x2);
        z2=(v1+v4-v2-v3)/ep3;
        for m=1:k
            s(m).hf(i,j)=z2(m);
            s(m).hf(j,i)=s(m).hf(i,j);
        end
        x1(j)=a(j);
        x2(j)=a(j);
        j=j+1;
    end
end
f=s; % return structure array containing all the hessian matrices of 
%      constraint functions
return;
end