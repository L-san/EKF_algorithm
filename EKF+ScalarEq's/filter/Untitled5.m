n = 100;
% Qq = [6.84210526315790e-10];
% Qw = [5.26315789473684e-11];
Qq = 9.737e-11;%0.555555555555556;
Qw = 2.395e-10;%8888888889;
Q = diag([Qq Qq Qq Qw Qw Qw]);
L = zeros(6,n);
% for i = 1:n
%     L(:,i) = main(strcat('data',int2str(i),'.mat'),Q);
% end
%Q = eye(6,6)*1e-4;
main("data_distb2.mat",Q,1);

