%%% Linearization around a given trimPoint
%% Choose trimpoint matrix, from simulink trim tool.

load('C:\Users\rapha\Desktop\MA2\Advanced control\surface\Workspace\AdvancedControlCourse\LQR_Raphael\OpTrim60.mat');

u0=[60 60]'; %To change accordingly to what's in OpTrimXX.Inputs
x0=OpTrim60.states.x;
% sys = linmod('OtterFullForLqr',OpTrim60.states.x,u0); % ('simulinkModelName','trimCond')

load('C:\Users\rapha\Desktop\MA2\Advanced control\surface\Workspace\AdvancedControlCourse\LQR_Raphael\LinSyst60.mat');

A = sys.a;
B = sys.b;
C = sys.c;
D = sys.d;

G = tf(ss(A,B,C,D));

%% LQR design.
Q = eye(size(sys.a,1));%Put higher to allow less error for the states. 0 for states in which we don't care about their error.
% 
Q(1,1)=1;
Q(6,6) = 1;%Yaw angle velocity

R = eye(size(sys.b,2));% Put higher to lower the use of the propellers
[F,P,CLP] = lqr(A,B,Q,R);
F(:,2:5)=0;
F(:,7:12)=0;
