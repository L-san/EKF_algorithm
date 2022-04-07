n = 100;
Qq = 4e-4;
Qw = 6e-4;
for i = 1:n
     main(strcat('data',int2str(i),'.mat'));
end

main("data_distb.mat",diag([Qq Qq Qq Qw Qw Qw]));

main("data_distb2.mat")