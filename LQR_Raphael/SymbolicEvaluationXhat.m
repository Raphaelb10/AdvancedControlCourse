function [A]=SymbolicEvaluationXhat(xhat_kplus1,Thrust)
%% Obtaining symbolic function
[F,X,U]=OtterFunction_EKF(); %xdot=f(x,u);
A=jacobian(F,X); %Compute A matrix of ekf, symbolic form


%% Eval xhat
% xhat_kplus1= OtterFunction_lqr(xhat_k,u_k,0,0,0,0);

A=subs(A, X, xhat_kplus1);
A=subs(A,U,Thrust);
A=double(A);% Cast it to double

end