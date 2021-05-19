%orbit_vec = [sma ecc inc raan aop ta w0]; %orbital elements vector
function attitudeVector = motionEquations(t,attitude)
global orbit_vec satellite mu
%quaternion and angular velocity in the body system
q = attitude(1:4);
w_b = attitude(5:7);
r = attitude(8:10);
V = attitude(11:13);

%mean motion
Aorb2b = quat2dcm(attitude(1:4)');
w0_b = Aorb2b*[0; orbit_vec.w0; 0]; 

%angular velocity in the orbital system
w_orb = w_b+w0_b;

%poisson's eq-s
q_dot(1,:) = -0.5*dot(q(2:4), w_orb);%scalar part
q_dot(2:4,:) = 0.5 * (q(1) * w_orb + cross(q(2:4), w_orb));%vector part

%external torques
M = calculateTorques(attitude);

%Eulers eq-s
w_dot = satellite.I \ (M - cross(w_orb, satellite.I * w_orb));

%Orbital motion
k = -mu/norm(r)^3;

attitudeVector = [q_dot; w_dot; V; k*r];
end