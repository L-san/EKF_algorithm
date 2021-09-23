function z = getHiddenState(attitude ,Borb)
mu = 398600e+9;
q(:,1) = attitude(1:4);
w(:,1) = attitude(5:7);
r(:,1) = attitude(8:10);
V(:,1) = attitude(11:13);
w0 = [0;sqrt(mu/(norm(r))^3);0];
B_Ax = Borb(1); B_Ay = Borb(2); B_Az = Borb(3);
q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);
Aorb2b = quat2DCM(q');
%B_b = quat_mult(quat_mult(q,[0;Borb]),quat_conj(q));

Bb1 = 0;
Bb2 = 2*(B_Az*q1*q3+B_Ay*q2*q3-B_Ay*q1*q4+B_Az*q2*q4)+B_Ax*(q1^2+q2^2-q3^2-q4^2);
Bb3 = 2*(-B_Az*q1*q2+B_Ax*q2*q3+B_Ax*q1*q4+B_Az*q3*q4)+B_Ay*(q1^2-q2^2+q3^2-q4^2);
Bb4 = 2*(B_Ay*q1*q2-B_Ax*q1*q3+B_Ax*q2*q4+B_Ay*q3*q4)+B_Az*(q1^2-q2^2-q3^2+q4^2);
B_b = [Bb1; Bb2; Bb3; Bb4];
%qdot = 0.5*quat_mult(q,w);
W = w+Aorb2b*w0;
z = [B_b; W];
end