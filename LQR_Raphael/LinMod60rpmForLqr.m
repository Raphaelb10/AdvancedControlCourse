%%% Linearization around a given trimPoint
%% Choose trimpoint matrix, from simulink trim tool.
addpath(genpath('C:\Users\rapha\Desktop\MA2\Advanced control\surface\MSS-master\MSS-master')); %%Add mss-master to path
addpath(genpath("C:\Users\rapha\Desktop\MA2\Advanced control\surface\Workspace\MSSfor2020"));
load('C:\Users\rapha\Desktop\MA2\Advanced control\surface\Workspace\AdvancedControlCourse\LQR_Raphael\OpTrim60.mat');

u0=[60 60]'; %To change accordingly to what's in OpTrimXX.Inputs
x0=OpTrim60.states.x;
% sys = linmod('OtterFullForLqr',OpTrim60.states.x,u0); % ('simulinkModelName','trimCond')

load('C:\Users\rapha\Desktop\MA2\Advanced control\surface\Workspace\AdvancedControlCourse\LQR_Raphael\LinSyst60.mat');

A = sys.a;
B = sys.b;
Cfull = sys.c;
Dfull = sys.d;

G = tf(ss(A,B,Cfull,Dfull));
%% Changing C matrix
C=[1 0 0 0 0 0 0 0 0 0 0 0;
   0 0 0 0 0 1 0 0 0 0 0 0];
D=zeros(2,2);
%% LQR design.
Q = eye(size(sys.a,1));%Put higher to allow less error for the states. 0 for states in which we don't care about their error.
% 
Q(1,1)=1;
Q(6,6)=1;%Yaw angle velocity

R = eye(size(sys.b,2));% Put higher to lower the use of the propellers
[F,P,CLP] = lqr(A,B,Q,R);
F(:,2:5)=0;
F(:,7:12)=0;
%% Outerloop control for position
sys=ss(A,B,C,D);
Dcgain=dcgain(sys);
N=inv(Dcgain);
%sI_tf=tf([1 0],1)*eye(12)
% N=(C*inv((eye(12)-(A-B*F)))*B)'%end up with 12x2 instead
% N(:,2:5)=0;
% N(:,7:12)=0