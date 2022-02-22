function x = getAttitudeVector(attitude, dt)
att_vec = motionEquations(attitude);
q_dot = att_vec(1:4,:);
w_dot = att_vec(5:7,:);

q = attitude(1:3,:);
w = attitude(4:6,:);

q = q+q_dot(2:4)*dt;
w = w+w_dot*dt;

x = [q; w];
end