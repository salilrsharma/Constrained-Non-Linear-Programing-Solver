function [f] = GE1(A,b)
% Objective: Implements gaussian elemination to solve Ax=b
%-----------------------------------------------------------------------
% [x] = GE1(A,b)
% where A=LHS of Ax=b 
%       b=RHS of Ax=b
%-----------------------------------------------------------------------
% Output: [f]= Column array of solution to a system of equations
%-----------------------------------------------------------------------

% Code by:
% Salil Sharma
% For the project implementation in IE 538 course
% Spring 2017
%-----------------------------------------------------------------------

N = length(b); %length of b vector
for c=1:(N-1) %loop through columns
    %swap rows 
    [gues,ind] = max(abs(A(c:end,c)));
    ind=ind+c-1;
    temp  = A(c,:);
    A(c,:) = A(ind,:);
    A(ind,:)  = temp;
    temp = b(c);
    b(c)= b(ind);
    b(ind) = temp;
    %all rows lying below the diagonal
    for r =(c+1):N %loop through rows
        d = A(r,c)/A(c,c); % compute d value
        A(r,c:end) = A(r,c:end)-d*A(c,c:end) ;
          b(r)          = b(r)-d*b(c);
    end
end 
for r=N:-1:1 %back substitution
x(r) = b(r);
    for i=(r+1):N
          x(r) = x(r)-A(r,i)*x(i);
    end
x(r) = x(r)/A(r,r); %solve for x(i)
end
f = x'; % returns solution [x]
return;