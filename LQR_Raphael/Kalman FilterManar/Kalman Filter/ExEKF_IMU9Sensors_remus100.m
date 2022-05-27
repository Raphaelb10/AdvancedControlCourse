% ExEKF Discrete-time extended Kalman filter (EKF) implementation demonstrating
% how the "predictor-corrector representation" can be applied to the
% nonlinear model: remus100

%% USER INPUTS
clear all
f_s    = 10;    % sampling frequency [Hz]
f_IMU = f_s;     % IMU Reading frequency   
f_gnss= 0.1; % GNSS position measurement frequency [Hz]

% simulation parameters
N  = 150;		  % no. of iterations
h  = 1/f_s; 	  % sampling time: h  = 1/f_s (s) 
h_gnss = 1/f_gnss;  

% remus model xdot, states, input= remus100Sym12States xdot=F(x,u) 
[F,X,U]=remus100Sym12States();% obtaining the function xdot= F(X,U) in symbolic variables 
A= jacobian(F,X);
% initial values for x and u
x = [1 1 1 1 1 1 1 1 1 1 1 1 ]';	% could be any value, choosing a zero vector creates problems with Jacobian
u = [0 0 0]';

% initialization of Kalman filter
x_prd = x;  % You know the initial state   
% P_prd is 12x12 
P_prd = zeros(12); % should be zero because you are sure of your initial state 
Qd = eye(12); % should be of square dimension as that of state: 12x12
Rd = eye(9); % should be of square dimension as that of output: 9x9  

% noise standard deviations
sigma_gyro=0.1; % rad/s
sigma_accel=0.1; % m/s/s
sigma_magn= 1*pi/180;% rad

% disturbance standard deviations 
sigma_angular_vel= 0.1; % rad/s
sigma_linear_vel= 0.1;% m/s
sigma_pos= 0.5;% m
sigma_euler_ang= 5*pi/180;% rad

% The IMU measures 9 outputs: linear accelerations, angular velocities and
% euler angles 

%% MAIN LOOP
simdata = zeros(N+1,25); % table of simulation data:1 for time, 12 for real states and 12 for estimated states
ydata= zeros(N+1,10) ; % table of measurement data: 1 for time, 9 for outputs 

for i=1:N+1
   t = (i-1) * h;                          % time (s)             

   % Plant( non-linear continuous model) 
   
   % input
   if 5*(i-1)<1500 % have not reached saturation? 
   u = [1.2 * sin(0.6*t) 0 5*(i-1)]'; 
   else
   u = [1.2 * sin(0.6*t) 0 1500]';    
   end
   % disturbance 
          w =[sigma_linear_vel * randn(1);
             sigma_linear_vel * randn(1);
             sigma_linear_vel * randn(1);
             sigma_angular_vel * randn(1);
             sigma_angular_vel * randn(1);
             sigma_angular_vel * randn(1);
             sigma_pos * randn(1);
             sigma_pos * randn(1);
             sigma_pos * randn(1);
             sigma_euler_ang * randn(1);
             sigma_euler_ang * randn(1);
             sigma_euler_ang * randn(1)];
         
   % w has the dimension of the state vector 
   [f, non_sense]= remus100(x,u); % xdot=f(x,u)+ w    
    x_dot = f + w;                      % dx/dt = f +  w  
   
  

   % IMU measurements at freq= sampling freq
   % The IMU measures 9 outputs: linear accelerations, angular velocities and
    % euler angles. The linear accelerations should be integrated to the
    % linear velocities. Note that these are along the BODY frame 
        
    
    if i==1
        
    % Initital outputs, note that the first three are not initial readings
    % of IMU beacause first three initial readings would be accelerations
    
    y =[ x(1) ;
         x(2) ;
         x(3) ;
         x(4) + sigma_gyro * randn(1);
         x(5) + sigma_gyro* randn(1);
         x(6) + sigma_gyro * randn(1);
         x(10) + sigma_magn * randn(1);
         x(11) + sigma_magn * randn(1);
         x(12) + sigma_magn * randn(1)];
     
    else
         y= [ y_next(1);
              y_next(2);
              y_next(3);
         x(4) + sigma_gyro * randn(1); % y(4:6) is gyroscope
         x(5) + sigma_gyro * randn(1);
         x(6) + sigma_gyro * randn(1);
         x(10) + sigma_magn * randn(1);% y(7:9) is from magnetometer
         x(11) + sigma_magn * randn(1);
         x(12) + sigma_magn * randn(1)];
        
        
    end
     % Euler integration (k+1)
        ax= x_dot(1) + sigma_accel * randn(1);
        ay= x_dot(2) + sigma_accel * randn(1);
        az= x_dot(3) + sigma_accel * randn(1);
    
        y_next(1)=y(1)+h*ax;
        y_next(2)=y(2)+h*ay;
        y_next(3)=y(3)+h*az;

     
   
    ydata(i,:)=[ t, y']; 
       % Cd is 9x12 
    Cd = [1 0 0 0 0 0 0 0 0 0 0 0;
          0 1 0 0 0 0 0 0 0 0 0 0;
          0 0 1 0 0 0 0 0 0 0 0 0;
          0 0 0 1 0 0 0 0 0 0 0 0;
          0 0 0 0 1 0 0 0 0 0 0 0;
          0 0 0 0 0 1 0 0 0 0 0 0;
          0 0 0 0 0 0 0 0 0 1 0 0;
          0 0 0 0 0 0 0 0 0 0 1 0;
          0 0 0 0 0 0 0 0 0 0 0 1];
 
    
   % KF is of dimesnions 12x9
   % KF gain      
   K = P_prd * Cd' * inv( Cd * P_prd * Cd' + Rd ); % P_prd is P_k+1 
        
   % corrector   
   % correction of P_prd
   IKC = eye(12) - K*Cd;
   P_hat = IKC * P_prd * IKC' + K * Rd * K';% In this step in the notes, we did not take the effect of R 
   
   % correction of state
   eps = y - Cd * x_prd; %x_prd is x_k+1 hat, Cd*x_prd is your open loop estimation of the states that are outputs
   x_hat = x_prd + K * eps; % corrected value of the estimated state   
   
   % store simulation data in a table   
   simdata(i,:) = [t x' x_hat'];    

   % discrete-time extended KF-model
   %f_hat is the evaluation of f at x_hat
   f_hat = remus100(x_hat,u);
   %f_d is the discrete form of f , evaluated at x_hat
   f_d   = x_hat + h * f_hat;
      
   % Predictor (k+1)  
   % Ad = I + h * A due to discretization of co
   % where A = df/dx is linearized about x = x_hat
   
    
   
    A=subs(A, X, x_hat);
    A=subs(A,U,u);
    A=double(A);% Ad is A_k in the notes 
    Ad= eye(12) + h*A;

   
   x_prd = f_d;
   P_prd = Ad * P_hat * Ad' + Qd;
   
   % Euler integration (k+1)
   x = x + h * x_dot; % continuous time states (real states)
   counter= i
end

%% PLOTS
close all
t     = simdata(1:i,1); 
surge     = simdata(1:i,2); 
surge_hat = simdata(1:i,14); 
x_pos    = simdata(1:i,8); 
x_pos_hat = simdata(1:i,20);

sway     = simdata(1:i,3); 
sway_hat = simdata(1:i,15); 
y_pos    = simdata(1:i,9); 
y_pos_hat = simdata(1:i,21);

heave     = simdata(1:i,4); 
heave_hat = simdata(1:i,16); 
z_pos    = simdata(1:i,10); 
z_pos_hat = simdata(1:i,22);

angular_p     = simdata(1:i,5); 
angular_p_hat = simdata(1:i,17); 
roll_angle    = simdata(1:i,11); 
roll_angle_hat = simdata(1:i,23);

angular_q     = simdata(1:i,6); 
angular_q_hat = simdata(1:i,18); 
pitch_angle    = simdata(1:i,12); 
pitch_angle_hat = simdata(1:i,24);

angular_r     = simdata(1:i,7); 
angular_r_hat = simdata(1:i,19); 
yaw_angle    = simdata(1:i,13); 
yaw_angle_hat = simdata(1:i,25);


t_m = ydata(:,1); % time 
surge_m = ydata(:,2); % state measurement 
sway_m = ydata(:,3); % state measurement 
heave_m = ydata(:,4); % state measurement 
angular_p_m = ydata(:,5); % state measurement 
angular_q_m = ydata(:,6); % state measurement 
angular_r_m = ydata(:,7); % state measurement 
roll_angle_m = ydata(:,8); % state measurement 
pitch_angle_m = ydata(:,9); % state measurement 
yaw_angle_m = ydata(:,10); % state measurement 



clf
subplot(611),plot(t_m,surge_m,'xb',t,surge_hat,'r')
xlabel('time (s)'),title('surge velocity'),grid
legend(['Measurement of surge at ', num2str(f_s), ' Hz'],...
    ['Estimate of u at ', num2str(f_s), ' Hz']);

subplot(612),plot(t,x_pos,'b',t,x_pos_hat,'r')
xlabel('time (s)'),title('x position'),grid
legend(['True x position at ', num2str(f_s), ' Hz'],...
    ['Estimate of x at ', num2str(f_s), ' Hz']);

subplot(613),plot(t_m,sway_m,'xb',t,sway_hat,'r')
xlabel('time (s)'),title('sway velocity'),grid
legend(['Measurement of sway at ', num2str(f_s), ' Hz'],...
    ['Estimate of sway at ', num2str(f_s), ' Hz']);

subplot(614),plot(t,y_pos,'b',t,y_pos_hat,'r')
xlabel('time (s)'),title('y position'),grid
legend(['True y position at ', num2str(f_s), ' Hz'],...
    ['Estimate of y at ', num2str(f_s), ' Hz']);

subplot(615),plot(t_m,heave_m,'xb',t,heave_hat,'r')
xlabel('time (s)'),title('heave velocity'),grid
legend(['Measurement of heave at ', num2str(f_s), ' Hz'],...
    ['Estimate of heave at ', num2str(f_s), ' Hz']);

subplot(616),plot(t,z_pos,'b',t,z_pos_hat,'r')
xlabel('time (s)'),title('z position'),grid
legend(['True z position at ', num2str(f_s), ' Hz'],...
    ['Estimate of z at ', num2str(f_s), ' Hz']);

figure()

subplot(611),plot(t_m,angular_p_m,'xb',t,angular_p_hat,'r')
xlabel('time (s)'),title('Angular velocity p'),grid
legend(['Measurement of p at ', num2str(f_s), ' Hz'],...
    ['Estimate of p at ', num2str(f_s), ' Hz']);

subplot(612),plot(t,roll_angle_m,'xb',t,roll_angle_hat,'r')
xlabel('time (s)'),title('Roll angle'),grid
legend(['Measurement of roll at ', num2str(f_s), ' Hz'],...
    ['Estimate of p at ', num2str(f_s), ' Hz']);

subplot(613),plot(t_m,angular_q_m,'xb',t,angular_q_hat,'r')
xlabel('time (s)'),title('Angular velocity q'),grid
legend(['Measurement of q at ', num2str(f_s), ' Hz'],...
    ['Estimate of q at ', num2str(f_s), ' Hz']);

subplot(614),plot(t,pitch_angle_m,'xb',t,pitch_angle_hat,'r')
xlabel('time (s)'),title('Pitch angle'),grid
legend(['Measurement of roll at ', num2str(f_s), ' Hz'],...
    ['Estimate of pitch angle at ', num2str(f_s), ' Hz']);

subplot(615),plot(t_m,angular_r_m,'xb',t,angular_r_hat,'r')
xlabel('time (s)'),title('Angular velocity p'),grid
legend(['Measurement of r at ', num2str(f_s), ' Hz'],...
    ['Estimate of r at ', num2str(f_s), ' Hz']);

subplot(616),plot(t,yaw_angle_m,'xb',t,yaw_angle_hat,'r')
xlabel('time (s)'),title('Yaw angle'),grid
legend(['Measurement of roll at ', num2str(f_s), ' Hz'],...
    ['Estimate of yaw angle at ', num2str(f_s), ' Hz']);

