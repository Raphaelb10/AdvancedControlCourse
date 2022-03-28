# Project on autonomous boat "Otter".
## How to test something.
The reference files are OtterFullModel and OtterFunction. 
To try something (e.g. linearize, simplify, control), make a copy then track it.

### TO DO list
-Something is wrong with the angle and the current. When we modify the current angle, the yaw angle seems to be not of the right sign.
Even with the non-modified model this happens. To reproduce, just try to make the boat turn by two means : by putting less thrust on one propeller, see the direction in which it'll turn. 
then try to put some current in the same direction. The XY plot will turn in the same direction, but if you look at the Euleur angle plot, you'll see the sign is not the same in the two situations.
Even thought the XY plot turns in the same direction. Any suggestions ?

Also it's strange that in fossen's model the XY were inverted in the XY plot. On documentation it is X at the top and Y at the bottom. I kept it as it is as everything seems normal.

-Linearize model with two angles set to 0. 
- Apply LQR on it


### File descriptions
-OtterFullModel & OtterFunction : Main reference for the full model. Do not modify
-OtterFixedAngleModel & Function : Simplified version where phi and theta are set to 0 (as eq) to simplify the dynamic.

Based on :
T. I. Fossen and T. Perez (2004). Marine Systems Simulator (MSS)
URL: https://github.com/cybergalactic/MSS