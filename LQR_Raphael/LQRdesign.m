%% LQR design.
Q = eye(size(sys.a,1));%Put higher to allow less error for the states. 0 for states in which we don't care about their error.
% 
Q(1,1)=10000;
Q(6,6)=10000;%Yaw angle velocity
Q(2,2)=0.0;
Q(3,3)=0.0;
Q(4,4)=0.0;
Q(9,9)=0.0;
Q(10,10)=0.0;
Q(11,11)=0.0;
Q(12,12)=0.0;

R = 10*eye(size(sys.b,2));% Put higher to lower the use of the propellers

[F,P,CLP] = lqr(A,B,Q,R);

F

F(:,2:5)=0;
F(:,7:12)=0;

F
