%%% Linearization around a given trimPoint
%% Choose trimpoint matrix, from simulink trim tool.

load('C:\Users\rapha\Desktop\MA2\Advanced control\surface\Workspace\AdvancedControlCourse\LQR_Raphael\OpTrim60.mat');

u0=[60 60]'; %To change accordingly to what's in OpTrimXX.Inputs
x0=OpTrim60.states.x;
sys = linmod('OtterFullForLqr',OpTrim60.states.x,u0); % ('simulinkModelName','trimCond')


A = sys.a;
B = sys.b;
C = sys.c;
D = sys.d;

G = tf(ss(A,B,C,D));

%% LQR design.

