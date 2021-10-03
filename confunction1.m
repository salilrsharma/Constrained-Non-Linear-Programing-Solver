function [c,ceq] = confunction1(x) % Nonlinear inequality constraints
%ceq = [x(1)^2+x(2)^2-4;(x(1)-1)^2+2*x(2)^2-4]; % Nonlinear equality constraints
%ceq=[([x(1),x(2),x(3),x(4)]*[x(1);x(2);x(3);x(4)])-4;([x(1)-1,x(2),x(3),x(4)]*[1,0,0,0;0,2,0,0;0,0,1,0;0,0,0,2]*[x(1)-1;x(2);x(3);x(4)])-4];
c=[];
ceq(1) = ([x(1),x(2),x(3)]*[x(1);x(2);x(3)])-4; % Nonlinear equality constraints
ceq(2)=([x(1)-1,x(2),x(3)]*[1,0,0;0,2,0;0,0,1]*[x(1)-1;x(2);x(3)])-4;
end