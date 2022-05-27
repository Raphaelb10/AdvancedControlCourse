function [xdot,x,u]= remus100Sym9States() 
% This function returns the symbolic form of xdot in terms of the state
% vector(9x1) and the inputs(3x1)


% [xdot,U] = remus100(x,ui,Vc,betaVc) returns the time derivative of the 
% state vector: x = [ u v w p q r x y z phi theta psi ]' and speed U in m/s  
% (optionally) for the Remus 100 autonomous underwater vehicle (AUV). The 
% length of the AUV is L = 1.7 m, while the state vector is defined as:
%
%  u:       surge velocity          (m/s)
%  v:       sway velocity           (m/s)
%  w:       heave velocity          (m/s)
%  p:       roll rate               (rad/s)
%  q:       pitch rate              (rad/s)
%  r:       yaw rate                (rad/s)
%  x:       North position          (m)
%  y:       East position           (m)
%  z:       downwards position      (m)
%  phi:     roll angle              (rad)
%  theta:   pitch angle             (rad)
%  psi:     yaw angle               (rad)
%
% The control inputs are one single rudder, two stern planes and one single 
% propeller:
%
%  ui = [ delta_r delta_s n ]'  where
%
%    delta_r:   rudder angle (rad)
%    delta_s:   stern plane angle (rad) 
%    n:         propeller revolution (rpm)
%
% The last arguments Vc and betaVc are optional arguments for ocean current 
% speed and direction expressed in NED.
%
% Author:    Thor I. Fossen
% Date:      27 May 2021
% Revisions: 24 Aug 2021, The ocean current is now expressed in NED 
%            10 Oct 2021, Increased the parasetic drag (alpha = 0) to 0.2
%            21 Oct 2021, implay61.m is called using the relative velocity
%            30 Dec 2021, Added the time derivative of the current velocity
%            01 Feb 2022, Updated lift and drag forces

% Check of input and state dimensions

% symbolic inputs
x = sym('x',[12,1],'real'); %[nu;eta]
u = sym('u',[3,1],'real'); % input 


% if (length(x) ~= 12),error('x-vector must have dimension 9!'); end
% if (length(ui) ~= 3),error('u-vector must have dimension 3!'); end


% Constants
mu = 63.446827;         % Lattitude for Trondheim (deg)
g_mu = gravity(mu);     % gravity vector (m/s2)
rho = 1026;             % density of water (m/s2)

% State vectors and control inputs
nu = x(1:6);
% /!\ 
eta = x(7:9);
% /!\ 
delta_r = u(1);        % tail rudder (rad)
delta_s = u(2);        % stern plane (rad)
n = u(3)/60;           % propeller revolution (rps)

% Ocean currents expressed in BODY
% /!\ 
    Vc = 0;
    betaVc = 0;
% /!\ 

% /!\ 
u_c = Vc * cos( betaVc - eta(3) );                               
v_c = Vc * sin( betaVc - eta(3) );  
% /!\ 

nu_c = [u_c v_c 0 0 0 0]';                  % ocean current velocities
Dnu_c = [nu(6)*v_c -nu(6)*u_c 0 0 0 0]';    % time derivative of nu_c

% Amplitude saturation of control signals
max_ui = [30*pi/180 70*pi/180  1500/60]';   % deg, deg, rps

% Relative velocities, speed and angle of attack
nu_r = nu - nu_c;                                 % relative velocity
alpha = atan2( nu_r(3), nu_r(1) );                % angle of attack (rad)
U_r = sqrt( nu_r(1)^2 + nu_r(2)^2 + nu_r(3)^2 );  % relative speed (m/s)
U  = sqrt( nu(1)^2 + nu(2)^2 + nu(3)^2 );         % speed (m/s)

% AUV model parameters (Section 8.4.2)
L_auv = 1.7;                        % AUV length (m)
D_auv = 0.19;                       % AUV diamater (m)
CD_0 = 0.2;                         % parasitic drag
S = 0.7 * L_auv * D_auv;            % S = 70% of rectangle L_auv * D_auv
a = L_auv/2;                        % semi-axes
b = D_auv/2;                  
r_bg = [ 0 0 0.02 ]';               % CG w.r.t. to the CO
r_bb = [ 0 0 0 ]';                  % CB w.r.t. to the CO

% Added moment of inertia in roll
r44 = 0.3;               % A44 = r44 * Ix

% Propeller data
D_prop = 0.055;          % propeller diameter
KT = 0.4739;             % [KT, KQ] = wageningen(0,1,0.8,3)
KQ = 0.0730;

% Tail rudder (single)
CL_delta_r = 0.5;        % rudder lift coefficient
A_r = 0.10 * 0.05;       % rudder area (m2)
x_r = -a;                % rudder x-position (m)

% Stern plane (double)
CL_delta_s = 0.7;        % stern-plane lift coefficient
A_s = 2 * 0.10 * 0.05;   % stern-plane area (m2)
x_s = -a;                % stern-plane z-position (m)

% Low-speed linear damping matrix parameters: D * exp(-3 * U_r)
T1 = 20;                 % time constant in surge (s)
T2 = 20;                 % time constant in sway (s)
zeta4 = 0.3;             % relative damping ratio in roll
zeta5 = 0.8;             % relative damping ratio in pitch
T6 = 5;                  % time constant in yaw (s)

% Rigid-body mass and hydrodynamic added mass
[MRB,CRB] = spheroid(a,b,nu(4:6),r_bg);
[MA,CA] = imlay61(a, b, nu_r, r44);

% Nonlinear quadratic velocity terms in pitch and yaw (Munk moments) 
% are set to zero since only linear damping is used
CA(5,1) = 0;   
CA(5,4) = 0;
CA(6,1) = 0;
CA(6,2) = 0;

M = MRB + MA;
C = CRB + CA;
m = MRB(1,1); W = m*g_mu; B = W;

% Dissipative forces and moments
%intermediate D matrix, because D matrix is initialized as a double matrix and it can't be converted to a symbolic matrix afterwards
D_int = Dmtrx([T1 T2 T6],[zeta4 zeta5],MRB,MA,[W r_bg' r_bb']);

% /!\ 
% D(1,1) = D(1,1) * exp(-3*U_r);   % vanish at high speed where crossflow
% D(2,2) = D(2,2) * exp(-3*U_r);   % drag and lift/drag dominates
% D(6,6) = D(6,6) * exp(-3*U_r);

D = [D_int(1,1)*exp(-3*U_r), D_int(1,2), D_int(1,3), D_int(1,4), D_int(1,5),D_int(1,6);
    D_int(2,1), D_int(2,2)*exp(-3*U_r), D_int(2,3), D_int(2,4), D_int(2,5),D_int(2,6);
    D_int(3,1), D_int(3,2), D_int(3,3), D_int(3,4), D_int(3,5),D_int(3,6);
    D_int(4,1), D_int(4,2), D_int(4,3), D_int(4,4), D_int(4,5),D_int(4,6);
    D_int(5,1), D_int(5,2), D_int(5,3), D_int(5,4), D_int(5,5),D_int(5,6);
    D_int(6,1), D_int(6,2), D_int(6,3), D_int(6,4), D_int(6,5),D_int(6,6)*exp(-3*U_r)];

% /!\ 
tau_liftdrag = forceLiftDrag(D_auv,S,CD_0,alpha,U_r);
tau_crossflow = crossFlowDrag(L_auv,D_auv,D_auv,nu_r);

% Restoring forces and moments
% /!\ 
g = gvect(W,B,eta(2),eta(1),r_bg,r_bb);
% /!\ 

% kinematics
% /!\ 
J = eulerang(eta(1),eta(2),eta(3));
% /!\ 

% amplitude saturation of the control signals
% /!\ 
% if (abs(delta_r) > max_ui(1)), delta_r = sign(delta_r) * max_ui(1); end
% if (abs(delta_s) > max_ui(2)), delta_s = sign(delta_s) * max_ui(2); end
% if (abs(n)       > max_ui(3)), n = sign(n) * max_ui(3); end
% /!\ 

% Control forces and moments
X_prop = rho * D_prop^4 * KT * abs(n) * n;  % propeller thrust 
K_prop = rho * D_prop^5 * KQ * abs(n) * n;  % propeller-induced roll moment

% Horizontal- and vertical-plane relative speed
U_rh = sqrt( nu_r(1)^2 + nu_r(2)^2 );  
U_rv = sqrt( nu_r(1)^2 + nu_r(3)^2 );  

% Rudder and stern-plane drag
X_r = -0.5 * rho * U_rh^2 * A_r * CL_delta_r * delta_r^2; 
X_s = -0.5 * rho * U_rv^2 * A_s * CL_delta_s * delta_s^2;

% Rudder sway force 
Y_r = -0.5 * rho * U_rh^2 * A_r * CL_delta_r * delta_r;

% Rtern-plane heave force
Z_s = -0.5 * rho * U_rv^2 * A_s * CL_delta_s * delta_s;

% Generalized force vector
% /!\ 
tau = sym('tau',[6,1],'real'); 
% /!\ 
tau(1) = X_prop + X_r + X_s;
tau(2) = Y_r;
tau(3) = Z_s;
tau(4) = K_prop;
tau(5) = x_s * Z_s;
tau(6) = x_r * Y_r;

% /!\ 
velocities= J*nu; 
% /!\ 
% State-space model

% /!\ 
% inv(M) works with symbolic, unlike " \ " operator 
xdot = [...
  Dnu_c + inv(M)*(tau + tau_liftdrag + tau_crossflow - C * nu_r - D * nu_r  - g)
  velocities(4:6) ];
% /!\ 

end 
