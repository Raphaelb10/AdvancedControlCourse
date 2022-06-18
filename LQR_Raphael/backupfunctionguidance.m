function [u_d,psi_d] = fcn(ref_x,ref_y,x,y,psi)
%% Initialize the ref point we want to track
Npoint=size(ref_x,1); %Number of waypoints
k=1;
if(k<=Npoint)
    x_ref=ref_x(1);
    y_ref=ref_y(1);
    %% Compute desired yaw angle
    psi_d = atan2(y_ref-y, x_ref-x);  % Calculate the angle between the two points
    %% Error calculation
    error_psi=psi_d-psi;
    error_dist=sqrt((y_ref-y)^2+(x_ref-x)^2);

    %% Compute surge velocity wished depending on the boat orientation.
    if (error_dist > 1 && error_psi<5*(pi/180)) %When far away from set point and angle is already close to the wished one
        u_d=error_dist;
    else
        u_d=error_dist/100;%Move slowly when close to the point or if the angle error is too big.
%         k=k+1; %Go to next waypoint, as close enought to this one
    end
else
    u_d=0;%At the end don't move and keep same angle
    psi_d=psi;
end


