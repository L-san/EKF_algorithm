function JH = observationMatrix(attitude, Borb)
global orbit_vec
w0 = orbit_vec.w0;
q0 = attitude(1); q1 = attitude(2); q2 = attitude(3); q3 = attitude(4);
BAx = Borb(1); BAy = Borb(2); BAz = Borb(3);

E11 = BAx*q0+BAz*q2-BAy*q3;
E12 = BAx*q1+BAy*q2+BAz*q3;
E13 = BAz*q0+BAy*q1-BAx*q2;
E14 = -BAy*q0+BAz*q1-BAx*q3;  

E21 = BAy*q0+-BAz*q1+BAx*q3;
E22 = -BAy*q1-BAz*q0+BAx*q2;
E23 = BAy*q2+BAx*q1+BAz*q3;
E24 = BAx*q0+BAz*q2-BAy*q3; 

E31 = BAz*q0+BAy*q1-BAx*q2;
E32 = -BAz*q1+BAy*q0+BAx*q3;
E33 = -BAz*q2-BAx*q0+BAy*q3;
E34 = BAx*q1+BAy*q2+BAz*q3; 


G11 = -q3*w0;
G12 = q2*w0;
G13 = q1*w0;
G14 = -q0*w0;

G21 = q0*w0;
G22 = -q1*w0;
G23 = q2*w0;
G24 = -q3*w0;

G31 = -q1*w0;
G32 = -q0*w0; 
G33 = q3*w0;
G34 = q2*w0;

E = [0 0 0 0;
E11 E12 E13 E14;
E21 E22 E23 E24;
E31 E32 E33 E34];

F = zeros(4,3);

G = [G11 G12 G13 G14;
G21 G22 G23 G24;
G31 G32 G33 G34];

H = eye(3);

JH = [E, F;
      G, H];

end
