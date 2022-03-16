function x = getAttitudeVector(attitude, dt)
global q0
att_vec = motionEquations(attitude);
q_dot = att_vec(1:4,:);
w_dot = att_vec(5:7,:);

q = attitude(1:3,:);
w = attitude(4:6,:);

q = q+q_dot(2:4)*dt;
w = w+w_dot*dt;
q0 = q0 +q_dot(1)*dt;
x = [q; w];
end