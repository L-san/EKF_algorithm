function x = getAttitudeVector(attitude, dt)
att_vec = motionEquations(attitude);
q_dot = att_vec(1:4,:);
w_dot = att_vec(5:7,:);

q = attitude(1:4,:);
w = attitude(5:7,:);

q = q+q_dot*dt;
w = w+w_dot*dt;

x = [q; w];
end