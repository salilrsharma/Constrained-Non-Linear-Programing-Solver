options = optimoptions('fmincon','Algorithm','sqp');
fun = @(x)-x(1)+x(2)^2-x(3)^3+x(4)^4-x(5)^5;
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
nonlcon = @confunction1;
x0 = [0.2,0.2,0.1,0.1];
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)