function H = observationMatrix(attitude, Borb)
q0 = attitude(1); q1 = attitude(2); q2 = attitude(3); q3 = attitude(4);
wx = attitude(5); wy = attitude(6); wz = attitude(7);
Borb_x = Borb(1); Borb_y = Borb(2); Borb_z = Borb(3);

E21=Borb_x*q0+Borb_z*q2-Borb_y*q3;
E22=Borb_x*q1+Borb_y*q2+Borb_z*q3;
E23=Borb_z*q0+Borb_y*q1-Borb_x*q2;
E24=-Borb_y*q0+Borb_z*q1-Borb_x*q3;

E31=Borb_y*q0-Borb_z*q1+Borb_x*q3;
E32=-Borb_z*q0-Borb_y*q1+Borb_x*q2;
E33=Borb_x*q1+Borb_y*q2+Borb_z*q3;
E34=Borb_x*q0+Borb_z*q2-Borb_y*q3;

E41=Borb_z*q0+Borb_y*q1-Borb_x*q2;
E42=Borb_y*q0-Borb_z*q1+Borb_x*q3;
E43=-Borb_x*q0-Borb_z*q2+Borb_y*q3;
E44=Borb_x*q1+Borb_y*q2+Borb_z*q3;


E = [0 0 0 0 ;
E21 E22 E23 E24;
E31 E32 E33 E34;
E41 E42 E43 E44];

F = zeros(4,3);
L = zeros(3,4);
M = diag([1,1,1]);

H = [E, F;
    L, M];

end

