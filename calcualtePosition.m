function [answer, time, A_I, A_b, q_dot, W, q, W0]= calcualtePosition(delta_T,simT,Rcoeff,Qcoeff)
%%
%globals---------------------------------------------------------------
simConstT = simT/delta_T;
H = 117; %height above sea level,m
theta = 53.241505*pi/180;
lambda = 50.221245*pi/180;
dt = 1/100;
%%
%arrays--------------------------------------------------------------------
time = zeros(1,simConstT);
x_hat = zeros(8,simConstT);
q_dot = zeros(4,simConstT);
%%
%initial conditions--------------------------------------------------------
A_I = calculateGravity(H,theta,lambda);
A_I = A_I/norm(A_I);
q_dot(:,1) = [0 0 0 0];
Pcov = eye(8);
I = diag([0.01 0.01 0.05]);
data = parseData();
A_b = data(1:4,:);
A_b(:,1) = A_b(:,1)/norm(A_b(:,1));
W = data(5:8,:);
q = data(9:12,:);
W0 = data(13:end,:);
x_hat(:,1) = [1 0 0 0 W(:,1)']';
%--------------------------------------------------------------------------

%simulation
for i = 2:simConstT
   A_b(:,i) = A_b(:,i)/norm(A_b(:,i));
   A_b_dot = (A_b(:,i)-A_b(:,i-1))/dt;
   aux = getAttitudeVector(x_hat(:,i-1),I, delta_T);
   f = aux(1:8);
   q_dot(:,i) = aux(9:12);
   h = getHiddenState(x_hat(:,i-1), A_I);
   JF = stateTransitionMatrix(x_hat(:,i-1), I, delta_T);
   ff = JF*x_hat(:,i-1);
   JH = observationMatrix(x_hat(:,i-1), A_I, f);
   gg = JH*x_hat(:,i-1);
   [x_hat(:,i), Pcov] = kalmanReal(Rcoeff,Qcoeff, [A_b(:,i);  A_b_dot], Pcov, JF, JH, h, f);
   time(i) = time(i-1)+delta_T;
end
answer = x_hat;
end