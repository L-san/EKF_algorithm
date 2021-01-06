%orbit_vec = [sma ecc inc raan aop ta w0]; %orbital elements vector
function attitudeVector = motionEquations(t,attitude)
global orbit_vec
%quaternion and angular velocity in the body system
q0 = attitude(1);
q1 = attitude(2);
q2 = attitude(3);
q3 = attitude(4);

w_x = attitude(5);
w_y = attitude(6);
w_z = attitude(7);

%mean motion
w0_b = [0; orbit_vec(7); 0];    
w0_b_x = w0_b(1);
w0_b_y = w0_b(2);
w0_b_z = w0_b(3);

%angular velocity in the orbital system
w_orb_x = w_x-w0_b_x;
w_orb_y = w_y-w0_b_y;
w_orb_z = w_z-w0_b_z;

%poisson's eq-s
q_dot_0 = -0.5*(q1*w_orb_x+q2*w_orb_y+q3*w_orb_z);
q_dot_1 =  0.5*(q0*w_orb_x+q2*w_orb_z-q3*w_orb_y);
q_dot_2 =  0.5*(q0*w_orb_y +q3*w_orb_x-q1*w_orb_z);
q_dot_3 =  0.5*(q0*w_orb_z +q1*w_orb_y-q2*w_orb_x);

M = calculateTorques();%external torques
M_x = M(1); M_y = M(2); M_z = M(3);
Ix = 0.01; Iy = 0.01; Iz = 0.02;
%Eulers eq-s
w_dot_x = (M_x-(Iz-Iy)*w_y*w_z)/Ix;
w_dot_y = (M_y-(Ix-Iz)*w_z*w_x)/Iy;
w_dot_z = (M_z-(Iy-Ix)*w_x*w_y)/Iz;

attitudeVector = [q_dot_0 q_dot_1 q_dot_2 q_dot_3 w_dot_x w_dot_y w_dot_z]';
end