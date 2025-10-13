function velocityVector = ArrivalVelocity(velocity,lightSpeed)
%ArrivalVelocity gets the velocity vector from a virtual particule at zero  to
%arrive to the destination
%   Input:
%       velocity: particle velocity
%   Output:
%       velocityVector: the direction of the virtual particule reaching the destination
A=velocity(1);
B=velocity(2);
C=velocity(3);
E=lightSpeed^2;

% Solution to the follwing set of equations: solve (A + x)*t=D; (B+y)*t=0; (C+z)*t=0; Power[x,2]+Power[y,2]+Power[z,2]=E for x,y,z,t
%x = sqrt(-B^2 - C^2 + E) and y = -B and z = -C and t = (D (A - sqrt(-B^2 - C^2 + E)))/(A^2 + B^2 + C^2 - E) and A^2 + B^2 + C^2!=E

velocityVector = [sqrt(E - B^2 - C^2), -B, -C];

end