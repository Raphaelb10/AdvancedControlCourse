% ExEKF Discrete-time extended Kalman filter (EKF) implementation demonstrating
% how the "predictor-corrector representation" can be applied to the
% nonlinear model: remus100

%% USER INPUTS
clear all
close all 
f_s    = 20;    % sampling frequency [Hz]
f_IMU = f_s;     % IMU Reading frequency   
f_gnss= 1; % GNSS position measurement frequency [Hz]                


% simulation parameters
N  = 200;		  % no. of iterations
h  = 1/f_s; 	  % sampling time: h  = 1/f_s (s) 
h_gnss = 1/f_gnss; 

% remus model xdot, states, input= remus100Sym12States xdot=F(x,u) 
[F,X,U]=remus100Sym12States();% obtaining the function xdot= F(X,U) in symbolic variables 
A= jacobian(F,X); % A matrix in symbolic form 

% initial values for x and u
x = [1 1 1 1 1 1 1 1 1 1 1 1 ]';% could be any value, choosing a zero vecor creates problems with Jacobian
x_ideal=x;
u = [0 0 0]';

% noise standard deviations
sigma_gyro=0.01; % rad/s
sigma_accel=0.01; % m/s/s
sigma_magn= 1*pi/180;% rad
sigma_GNSS= 0.01; % m

% disturbance standard deviations( should not be much because I have a good
% model
sigma_angular_vel= 0.005; % rad/s/s
sigma_linear_vel= 0.005;% m/s/s
sigma_pos= 0.01;% m/s
sigma_euler_ang= 0.5*pi/180;% rad/s

%% initialization of Kalman filter
x_prd = x;  % You know the initial state 
x_hat = x;
% P_prd is 12x12 
P_hat = zeros(12); % should be zero because you are sure of your initial state 
Qd = 0.01*eye(12); % should be of square dimension as that of state: 12x12
Rd = 0.1*eye(9); % should be of square dimension as that of output: 9x9

% The IMU measures 9 outputs: linear accelerations, angular velocities and
% euler angles 
y= x(1:3); % Initialization of velocities for IMU 

%% MAIN LOOP
simdata = zeros(N+1,25); % table of simulation data:1 for time, 12 for real states and 12 for estimated states
ydata= zeros(N,10) ; % table of measurement data: 1 for time, 9 for outputs 
simdata(1,:) = [0 x' x_hat'];

for i=1:N+1
   %% time
    t = (i-1) * h;                          % time (s)             

   % Plant( non-linear continuous model) 
   
   %% input
   if 5*(i-1)<1500 % have not reached saturation? 
   u = [1.2 * sin(0.6*t) 0 5*(i-1)]'; 
   else
   u = [1.2 * sin(0.6*t) 0 1500]';    
   end
   
   %% disturbance 
   % w has the dimension of the state vector 
     w =[sigma_linear_vel * randn(3,1);
         sigma_angular_vel * randn(3,1);
         sigma_pos * randn(3,1);
         sigma_euler_ang * randn(3,1)];
             
   
   %% System dynamics
   [f, non_sense]= remus100(x,u); 
   [f_ideal, non_sense]= remus100(x_ideal,u);  
   x_dot = f+w;
   x_dot_ideal = f_ideal;
   
   % IMU measurements at freq= sampling freq
   % The IMU measures 9 outputs: linear accelerations, angular velocities and
   % euler angles. The linear accelerations should be integrated to the
   % linear velocities. Note that these are along the BODY frame
   

   
   %% Outputs    
   % Initital outputs, note that the first three are not initial readings
   % of IMU beacause first three initial readings would be accelerations
    
  
   % Euler integration (k+1)
   Acceleration = x_dot(1:3) + sigma_accel*randn(3,1);

   y= [y(1:3) + h*Acceleration;
       x(4:6) + sigma_gyro*randn(3,1);% y(4:6) is gyroscope
       x(10:12) + sigma_gyro*randn(3,1);...% y(7:9) is from magnetometer
       ];
   
   ydata(i,:)=[t, y']; 

%    % GNSS measurements are Z times slower than the sampling time
    if mod( t, h_gnss ) == 0
       y1 = x(7) + sigma_GNSS * randn(1); %x position
       y2 = x(8) + sigma_GNSS * randn(1); %y posiiton
       y3 = x(9) + sigma_GNSS * randn(1); %z posiiton 
      
       y= [y(1:6); y1; y2; y3;y(7:9)]; 
            
             
      Rd =0.1*eye(12);
      Sensors = [1,2,3,4,5,6,7,8,9,10,11,12]; %indices of the sensors with respect to their place in the state vector       
      Cd = zeros(12,12);
        
                  
    else
      Rd=0.1*eye(9);
      Sensors = [1,2,3,4,5,6,10,11,12]; %indices of the sensors with respect to their place in the state vector       
      Cd = zeros(9,12);
        
    end
    
for index = 1:length(Sensors)
    Cd(index,Sensors(index)) = 1;
end
    
   %% Prediction
   % discrete-time extended KF-model
   
   f_hat = remus100(x_hat,u); %f_hat is the evaluation of f at x_hat
   f_d   = x_hat + h * f_hat; %f_d is the discrete form of f , evaluated at x_hat
      
   % Predictor (k+1)  
   % Ad = I + h * A due to discretization of co
   % where A = df/dx is linearized about x = x_hat
   
   A=subs(A, X, x_hat);
   A=subs(A,U,u);
   A=double(A);% Ad is A_k in the notes 
   Ad= eye(12) + h*A;

   x_prd = f_d;
   P_prd = Ad * P_hat * Ad' + Qd;
  
   %% Correction of the prediction
   % KF is of dimesnions 12x9
   % KF gain      
   K = P_prd * Cd' /( Cd * P_prd * Cd' + Rd ); % P_prd is P_k+1 
        
   % corrector   
   
   % correction of P_prd
   IKC = eye(12) - K*Cd;
   P_hat = IKC * P_prd * IKC' + K * Rd * K';% In this step in the notes, we did not take the effect of R 
   
   % correction of state
   eps = y - Cd * x_prd; %x_prd is x_k+1 hat, Cd*x_prd is your open loop estimation of the states that are outputs
   x_hat = x_prd + K * eps; % corrected value of the estimated state   
   
   %% store values
   % store simulation data in a table   
   simdata(i,:) = [t x_ideal' x_hat']; 
   
   %% Integrate
   % Euler integration (k+1)
   x = x + h * x_dot;
   x_ideal = x_ideal + h * x_dot_ideal;
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




subplot(311),plot(t,x_pos,'b',t,x_pos_hat,'r')
xlabel('time (s)'),title('x position'),grid
legend(['True x position at ', num2str(f_s), ' Hz'],...
    ['Estimate of x at ', num2str(f_s), ' Hz']);

subplot(312),plot(t,y_pos,'b',t,y_pos_hat,'r')
xlabel('time (s)'),title('y position'),grid
legend(['True y position at ', num2str(f_s), ' Hz'],...
    ['Estimate of y at ', num2str(f_s), ' Hz']);

subplot(313),plot(t,z_pos,'b',t,z_pos_hat,'r')
xlabel('time (s)'),title('z position'),grid
legend(['True z position at ', num2str(f_s), ' Hz'],...
    ['Estimate of z at ', num2str(f_s), ' Hz']);

figure()
subplot(311),plot(t_m,surge_m,'xb',t,surge_hat,'r')
xlabel('time (s)'),title('surge velocity'),grid
legend(['Measurement of surge at ', num2str(f_s), ' Hz'],...
    ['Estimate of u at ', num2str(f_s), ' Hz']);



subplot(312),plot(t_m,sway_m,'xb',t,sway_hat,'r')
xlabel('time (s)'),title('sway velocity'),grid
legend(['Measurement of sway at ', num2str(f_s), ' Hz'],...
    ['Estimate of sway at ', num2str(f_s), ' Hz']);



subplot(313),plot(t_m,heave_m,'xb',t,heave_hat,'r')
xlabel('time (s)'),title('heave velocity'),grid
legend(['Measurement of heave at ', num2str(f_s), ' Hz'],...
    ['Estimate of heave at ', num2str(f_s), ' Hz']);



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

