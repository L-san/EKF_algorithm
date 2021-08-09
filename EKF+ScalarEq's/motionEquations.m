%orbit_vec = [sma ecc inc raan aop ta w0]; %orbital elements vector
function attitudeVector = motionEquations(t,attitude)
global orbit_vec satellite mu

%quaternion and angular velocity in the body system
q = attitude(1:4);
w_b = attitude(5:7);
r = attitude(8:10);
V = attitude(11:13);
w0 = orbit_vec.w0;
%dcm
Aorb2b = quat2DCM(attitude(1:4)');
Aorb2b12 = Aorb2b(1,2);
Aorb2b22 = Aorb2b(2,2);
Aorb2b32 = Aorb2b(3,2);
 
Ix = satellite.I(1,1); Iy = satellite.I(2,2); Iz = satellite.I(3,3);
w_bx = w_b(1); w_by = w_b(2); w_bz = w_b(3);
q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);

%mean motion
w0_bx = w0*Aorb2b12; w0_by = w0*Aorb2b22; w0_bz = w0*Aorb2b32;

%angular velocity in the orbital system
w_orb_x = w_bx-w0_bx;
w_orb_y = w_by-w0_by;
w_orb_z = w_bz-w0_bz;

%poisson's eq-s
q_dot_0 = -0.5*(q1*w_orb_x+q2*w_orb_y+q3*w_orb_z);
q_dot_1 = 0.5*(q0*w_orb_x+q2*w_orb_z-q3*w_orb_y);
q_dot_2 = 0.5*(q0*w_orb_y+q3*w_orb_x-q1*w_orb_z);
q_dot_3 = 0.5*(q0*w_orb_z+q1*w_orb_y-q2*w_orb_x);
q_dot = [q_dot_0; q_dot_1; q_dot_2; q_dot_3];

%external torques
M = calculateTorques(attitude);
Mx = M(1); My = M(2); Mz = M(3);

%Eulers eq-s
w_dot_x = Mx/Ix + (Iy - Iz)/Ix * w_by*w_bz;
w_dot_y = My/Iy + (Iz - Ix)/Iy * w_bx*w_bz;
w_dot_z = Mz/Iz + (Ix - Iy)/Iz * w_bx*w_by;
w_dot = [w_dot_x; w_dot_y; w_dot_z];


%Orbital motion
k = -mu/norm(r)^3;

attitudeVector = [q_dot; w_dot; V; k*r];
end