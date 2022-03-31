% x0 = [0 0 0 0 0 0 0 0 0 0 0 0]';
% u0 = [0 0]';
% y0 = [0 0]';
% [x_eq, u_eq, y_eq, dx_eq, options] = trim('otterSimulink_own',[],[],y0);
% sys = linmod('otterSimulink_own',x_eq,u_eq)


%clear all;
load('C:\Users\maxim\Desktop\cours\MA2\advanced theory\perso\trim_ship03.mat');
x0 = [0 0 7 0 0 0]'; %x, y, xdot, ydot, psi, psidot
u0 = [100 100]';
y0 = [0 0 7 0 0 0]';
ix = []; %Which state is fixed
iu = []; %Select which input is fixed
iy = [3]; %Select third output to be fixed
%dx0 = [7 0 0 0 0 0]';
%idx = [1;3];
u = [73.3303 73.3303]';
%[x_eq, u_eq, y_eq, dx_eq, options] = trim('Ship01',x0,[],y0,[],[],iy);
u=[12.1104 12.1104]';
%sys = linmod('ship03',op_ship03.states.x,u)
sys = linmod('ship04',trim_ship04.states.x,u)
%%
A = sys.a;
B = sys.b;
C = sys.c;
D = sys.d;

G = tf(ss(A,B,C,D));

CI = [0 0 0.5 0 0 0]';

%%
x0 = [0 0 0 0.5 1.57 0]'; %x, y, xdot, ydot, psi, psidot
%[x_eq, u_eq, y_eq, dx_eq, options] = trim('Ship01',x0,[],y0,[],[],iy);
u=[12.1104 12.1104]';
%sys = linmod('ship03',op_ship03.states.x,u)
sys = linmod('ship03',op_trim10.states.x,u)

CI = [0 0 0.5 0.5 pi/4 0]';

%%
x0 = [0 0 0 0.5 1.57 0]'; %x, y, xdot, ydot, psi, psidot
%[x_eq, u_eq, y_eq, dx_eq, options] = trim('Ship01',x0,[],y0,[],[],iy);
u=[0 0]';
%sys = linmod('ship03',op_ship03.states.x,u)
sys = linmod('ship04',op_trim1.states.x,u)

CI = [10 10 0 0 pi/4 0]';
%CI = op_trim8.states.x;
%%
%x0 = [5 0 0 0 0 0]'; %x, y, xdot, ydot, psi, psidot
%[x_eq, u_eq, y_eq, dx_eq, options] = trim('Ship01',x0,[],y0,[],[],iy);
u=[100 100]';
CI = [5 0 -1.03 -7.27e-24 -10 0 0.101 0.123 0.144 -1.94e-25 0.184 0]';
%sys = linmod('ship03',op_ship03.states.x,u)
sys = linmod('otter_6dof',CI,u)


%CI = op_trim8.states.x;
%%
%Q = [0 0 0
%LQR
Q = eye(size(sys.a,1));
Q(3:5,3:5) = 0; Q(9:11,9:11) = 0;
Q(6,6) = 1000;
R = eye(size(sys.b,2));
[F,P,CLP] = lqr(sys.a,sys.b,Q,R);
F3 = F;
F3(:,3:5) = 0; F3(:,9:11)=0;
%%
%Reachability
Reach = [B A*B A^2*B A^3*B A^4*B A^5*B]