clear variables; clc;
n = 10;
a = 1e-12; b = 1e-13;
Qw = a:(b-a)/(n-1):b; 
Qq = Qw;
dzeta = ones(n,n);
F.t = ones(n,n);
F.wx = ones(n,n);
F.wy = ones(n,n);
F.wz = ones(n,n);

for i = 1:n
    for j = 1:n
        Q = diag([Qw(i) Qw(i) Qw(i) Qq(j) Qq(j) Qq(j)]);
        [dzeta(i,j),f] = main("data_distb2.mat",Q,0);
        if(isnan(dzeta(i,j)))
            dzeta(i,j) = 0;
        end
        if(isnan(f))
            F.t(i,j) = 0.5;
            F.wx(i,j) = 2;
            F.wy(i,j) = 2;
            F.wz(i,j) = 2;
        else
            F.t(i,j) = f(1);
            F.wx(i,j) = f(2);
            F.wy(i,j) = f(3);
            F.wz(i,j) = f(4);
        end
    end
end
fig1 = figure;
imagesc(Qw,Qq,F.t);
contourf(Qw,Qq,F.t)
xlabel("Qw");
ylabel("Qq");
colormap(fig1,flipud(colormap(winter)))
colorbar;
title("Ф, градусы");
fig2 = figure; 
imagesc(Qw,Qq,dzeta);
contourf(Qw,Qq,dzeta)
xlabel("Qw");
ylabel("Qq");
colormap(fig2,winter);
colorbar;
title("n, %");

% figure;
% subplot(2,2,1); 
% imagesc(Qw,Qq,F.t);
% contourf(Qw,Qq,F.t)
% xlabel("Qw");
% ylabel("Qq");
% colormap winter;
% colorbar;
% title("Ф, градусы");
% subplot(2,2,2); 
% imagesc(Qw,Qq,F.wx);
% contourf(Qw,Qq,F.wx)
% xlabel("Qw");
% ylabel("Qq");
% colormap winter;
% colorbar;
% title("w_x, градусы/c");
% 
% subplot(2,2,3);
% imagesc(Qw,Qq,F.wy);
% contourf(Qw,Qq,F.wy)
% xlabel("Qw");
% ylabel("Qq");
% colormap winter;
% colorbar;
% title("w_y, градусы/c");
% 
% subplot(2,2,4);
% imagesc(Qw,Qq,F.wz);
% contourf(Qw,Qq,F.wz)
% xlabel("Qw");
% ylabel("Qq");
% colormap winter;
% colorbar;
% title("w_z, градусы/c");