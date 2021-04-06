function z = getHiddenState(attitude)
q = attitude(1:4);
w = attitude(5:8);

B_orb = Ain2orb

ax = A_I(2); ay = A_I(3); az = A_I(4);

A_b = quat_mult(quat_mult(q,[0 ax ay az]'),quatConj(q));

qdot = 0.5*quat_mult(q,w);
A_b_dot = quat_mult(quat_mult(qdot,[0 ax ay az]'),quatConj(q))+quat_mult(quat_mult(q,[0 ax ay az]'),quatConj(qdot));

z = [A_b' A_b_dot']';

end