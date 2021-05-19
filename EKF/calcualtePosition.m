function [answer, time]= calcualtePosition(dt,time,Rcoeff,Qcoeff,y)
%%
%globals---------------------------------------------------------------
N = length(time);
%%
%arrays--------------------------------------------------------------------
time = zeros(1,N);
x_hat = zeros(7,N);
%%
%initial conditions--------------------------------------------------------
Pcov = eye(7);
I = diag([0.01 0.01 0.05]);
x_hat(:,1) = y(1:7,1);
%additional noise----------------------------------------------------------
Aorb2b = quat2dcm(y(1:4,:)');
B_I = magneticField(y(8:10,:));
for j = 1:N
    V = y(11:13,j);
    r = y(8:10,j);
    Ain2orb = inv([V/norm(V), cross(r,V)/norm(cross(r,V)), r/norm(r)]);
    B_orb = Ain2orb*B_I(:,j);
    B_b(:,j) = Aorb2b(:,:,j)*B_orb;
end

B_b = awgn(B_b, 40,'measured');
w = awgn(y(5:7,:), 40,'measured');

%--------------------------------------------------------------------------
%simulation
for i = 2:N
   f = getAttitudeVector([x_hat(:,i-1);y(8:10,j);y(11:13,j)],dt);
   [h,B_orb] = getHiddenState([x_hat(:,i-1);y(8:10,j);y(11:13,j)], B_I(:,i));
   JF = stateTransitionMatrix(x_hat(:,i-1),dt);
  % ff = JF*x_hat(:,i-1);
   JH = observationMatrix(x_hat(:,i-1), B_orb);
  % gg = JH*x_hat(:,i-1);
   [x_hat(:,i), Pcov] = kalmanReal(Rcoeff,Qcoeff, [0; B_b(:,i);  w(:,i)], Pcov, JF, JH, h, f);
   time(i) = time(i-1)+dt;
end
answer = x_hat;
end