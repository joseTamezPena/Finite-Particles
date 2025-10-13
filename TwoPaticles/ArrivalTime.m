function [timetoArrival,velocityVector] = ArrivalTime(position,velocity,distance,lightSpeed)
%ArrivalTime gets the time from a virtual particule at position and velocity to
%arrive to the destination
%   Input:
%       position: particle position
%       velocity: particle velocity
%       distance: distance to destination
%   Output:
%       timetoArrival: the time the light particule will reach destination
%       velocityVector: the direction of the virtual particule reaching the destination
Avelocity = [0,0];
A=position(1);
C=position(2);
B=velocity(1);
D=velocity(2);
G=distance;
E=lightSpeed^2;

% Solution to the follwing set of equations: solve A+(B+x)*t=G;C+(D+y)*t=0;x^2+y^2=E for x,y,t
%x = (sqrt(-(A - G)^2 (A^2 D^2 - A^2 E - 2 A B C D - 2 A D^2 G + 2 A E G + B^2 C^2 + 2 B C D G - C^2 E + D^2 G^2 - E G^2)) + A C D - B C^2 - C D G)/(A^2 - 2 A G + C^2 + G^2)
%y = (A^3 (-D) + C sqrt(-(A - G)^2 (A^2 D^2 - A^2 E - 2 A B C D - 2 A D^2 G + 2 A E G + B^2 C^2 + 2 B C D G - C^2 E + D^2 G^2 - E G^2)) + A^2 B C + 3 A^2 D G - 2 A B C G - 3 A D G^2 + B C G^2 + D G^3)/((A - G) (A^2 - 2 A G + C^2 + G^2))
%t = -((A - G) (A^2 - 2 A G + C^2 + G^2))/(sqrt(-(A - G)^2 (A^2 D^2 - A^2 E - 2 A B C D - 2 A D^2 G + 2 A E G + B^2 C^2 + 2 B C D G - C^2 E + D^2 G^2 - E G^2)) + A^2 B - 2 A B G + A C D + B G^2 - C D G) 


%x = (-sqrt(-(A - G)^2 (A^2 D^2 - A^2 E - 2 A B C D - 2 A D^2 G + 2 A E G + B^2 C^2 + 2 B C D G - C^2 E + D^2 G^2 - E G^2)) + A C D - B C^2 - C D G)/(A^2 - 2 A G + C^2 + G^2) and 
%y = (A^3 (-D) - C sqrt(-(A - G)^2 (A^2 D^2 - A^2 E - 2 A B C D - 2 A D^2 G + 2 A E G + B^2 C^2 + 2 B C D G - C^2 E + D^2 G^2 - E G^2)) + A^2 B C + 3 A^2 D G - 2 A B C G - 3 A D G^2 + B C G^2 + D G^3)/((A - G) (A^2 - 2 A G + C^2 + G^2)) and 
%t = -((A - G) (A^2 - 2 A G + C^2 + G^2))/(-sqrt(-(A - G)^2 (A^2 D^2 - A^2 E - 2 A B C D - 2 A D^2 G + 2 A E G + B^2 C^2 + 2 B C D G - C^2 E + D^2 G^2 - E G^2)) + A^2 B - 2 A B G + A C D + B G^2 - C D G) and A^2 - 2 A G + C^2 + G^2!=0 and A - G!=0 and -sqrt(-(A - G)^2 (A^2 D^2 - A^2 E - 2 A B C D - 2 A D^2 G + 2 A E G + B^2 C^2 + 2 B C D G - C^2 E + D^2 G^2 - E G^2)) + A^2 B - 2 A B G + A C D + B G^2 - C D G!=0

dsqa = sqrt(-(A-G)^2*((A*D)^2-E*(A^2)-2*A*B*C*D-2*A*G*(D^2)+2*A*E*G+(B*C)^2+2*B*C*D*G-E*(C^2)+(D*G)^2-E*(G^2)));
denom = A^2+C^2+G^2-2*A*G;
dbias = A*C*D - B*(C^2) - C*D*G;

Avelocity(1) = (dsqa + dbias)/denom;

atime = (G-A)/(B + Avelocity(1));
if (atime < 0) 
    Avelocity(1) = (dbias - dsqa)/denom;
    atime = (G-A)/(B + Avelocity(1));
end
Avelocity(2) = -C/atime - D;
velm = sqrt(dot(Avelocity,Avelocity))/lightSpeed;
if (abs(velm-1.0)>1.0e-5)
    Avelocity(2) = sqrt(E - Avelocity(1)^2);
    atime = -C/(D + Avelocity(2));
    if (atime < 0)
        Avelocity(2) = -Avelocity(2);
        atime = -C/(D + Avelocity(2));
    end
end

timetoArrival = atime;
velocityVector = Avelocity;
end