function [optsol,optofv]=cnlp(x,a1,b1) 
% Objective: Solves non-linear optimization with equality constraints
%            using Sequential Quadratic Programing
%-----------------------------------------------------------------------
% [optsol,optofv]=cnlp(x,a1,b1)
% where x=Initial guess (column vector)
%       a1=objective function in the form of @(x)function
%       b1=Equality constraints separated by ';' in nx1 matrix
%-----------------------------------------------------------------------
% Output: optsol= optimal solution
%         optofv=Function value at the optimal point
%-----------------------------------------------------------------------

% Code by:
% Salil Sharma
% For the project implementation in IE 538 course
% Spring 2017
%-----------------------------------------------------------------------

xk=x; %Initial guess
test_function=a1; %test_function refers to the objective function
confunction=b1; % confunction refers to the equality constraints
xlength=length(x); % returns length of the initial guess
llengt=size(confunction(xk));
llength=llengt(1);% % returns length of the initial lambda
lk=ones(1,llength)'; % Initial lambda
tol=0.0000001; % Tolerance
err=1;
while err>tol % Loop terminates when error is within the specified tolerance
    %Break the lagrange into two: compute hessian separatley for objective
    %and constraints and at the end combine all the values
    M11=solveHessian(xk,test_function); % Hessian of the objective function at iterative opitn xk.
    w=LagHessian(xk,confunction); % Hessian of the constraints at xk
    w2=w(1).hf; %
    w3=size(w2);
    M13=zeros(w3(1));
    for n=1:llength
        M13=-lk(n)*w(n).hf+M13; % add all the constraints hessian 
    end
    M1=M11+M13; % M1 matrix comprises the hessian of the lagrange function
    M2=-LagGrad(xk,confunction); % gradient of the constraints in M2
    M3=-M2'; % M3 is transpose of -M2
    M4=zeros(llength); % create a matix of zero elements as M4
    M5=horzcat(M1,M2); %concatenate M1 and M2 
    M6=horzcat(M3,M4); %concatenate M3 and M4
    M=vertcat(M5,M6); % M matrix as in equation M.x=b
    %Create b matrix
    b1=-gradientfd(xk,test_function)'; % gradient of obj function at xk
    b2=-confunction(xk); % gradient of constraints at xk
    b=vertcat(b1,b2); %Concatenate b1 and b2
    r=GE1(M,b); % calls Gaussian elimination to solve Mx=b
    C=mat2cell(r,[xlength,llength]); % x matrix thus solved will contain xk value at the top and lambda updates at the bottom
    pk=cell2mat(C(1,1)); % retrieve xk values from C array
    xk=xk+pk; % Update xk
    lk=cell2mat(C(2,1)); % Retrieve lambda updates
    for p=1:llength % compute gradient of lagrange
        nm(:,p)=lk(p)*M2(:,p);
    end
    nm=sum(nm,2);
    err=norm(-b1+nm); %norm of lagrange gradient
end
optsol=xk; % Optimal solution
optofv=a1(xk); % Optimal function value
return;
end