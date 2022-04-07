function y = getAttitudeVector(attitude)
global q0
% k1 = motionEquations(attitude);
% k2 = motionEquations(attitude+[k1/2; zeros(6,1)]);
% k3 = motionEquations(attitude+[k2/2; zeros(6,1)]);
% k4 = motionEquations(attitude+[k3; zeros(6,1)]);

%y = attitude(1:7,:)+1/6*(k1+2*k2+2*k3+k4);
y1 =  [q0; attitude(1:6,:)] + motionEquations([q0; attitude])*0.01;

% if(norm(y1(2:4,:))>1)
%     q0 = 0;
% else
%     q0 = sqrt(1-norm(y1(2:4,:))^2);
% end

q0 = y1(1);
y = y1(2:7,:);
end