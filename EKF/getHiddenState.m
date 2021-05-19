function [z,B_orb] = getHiddenState(attitude,B_I)
mu = 398600e+9;
q(:,1) = attitude(1:4);
w(:,1) = attitude(5:7);
r(:,1) = attitude(8:10);
V(:,1) = attitude(11:13);
w0 = [0;sqrt(mu/(norm(r))^3);0];

Aorb2b = quat2dcm(q');
Ain2orb = inv([V/norm(V), cross(r,V)/norm(cross(r,V)), r/norm(r)]);
B_orb = Ain2orb*B_I;

B_b = quat_mult(quat_mult(q,[0;B_orb]),quat_conj(q));

%qdot = 0.5*quat_mult(q,w);
W = w+Aorb2b*w0;

z = [B_b; W];

end