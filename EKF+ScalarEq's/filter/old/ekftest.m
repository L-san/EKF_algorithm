clear variables; clc;
mu = 1;
n = 20;
precision = [1e-6 1e-6 1e-6 1e-7 1e-7 1e-7]';
Rw = 0:(1e20-1e18)/(n-1):1e20; 
Rq = 0:(1e20-1e18)/(n-1):1e20;
Q = zeros(6);

F.q1 = ones(n,n);
F.q2 = ones(n,n);
F.q3 = ones(n,n);
F.wx = ones(n,n);
F.wy = ones(n,n);
F.wz = ones(n,n);

sigma_b = [0.000131769397996411,0.000130058124526070,0.000570616647019740];
sigma_wb = [0.000192198568219992,0.000201278499705392,0.000216338468205052];
sigma_q = [5.45113088737962e-09,5.77278048768739e-09,5.62983975157216e-09];

Q = zeros(6,6);
P = diag([sigma_q, sigma_wb]);

for i = 1:n
    for j = 1:n
        R = diag([Rq(i) Rq(i) Rq(i) Rq(i) Rq(i) Rq(i)]);
        Q = diag([Rw(j) Rw(j) Rw(j) Rw(j) Rw(j) Rw(j)]);
        f = main(P,R,Q);
        if(isnan(f))
            F.q1(i,j) = 0;
            F.q2(i,j) = 0;
            F.q3(i,j) = 0;
            F.wx(i,j) = 0;
            F.wy(i,j) = 0;
            F.wz(i,j) = 0;
        else
            F.q1(i,j) = f(1);
            F.q2(i,j) = f(2);
            F.q3(i,j) = f(3);
            F.wx(i,j) = f(4);
            F.wy(i,j) = f(5);
            F.wz(i,j) = f(6);
        end
    end
end

figure;
subplot(2,3,1);
imagesc(Rw,Rq,F.q1);
%contourf(Rw,Rq,F.q1)
xlabel("Rw");
ylabel("Rq");
colormap winter;
colorbar;

subplot(2,3,2);
imagesc(Rw,Rq,F.q2);
%contourf(Rw,Rq,F.q2)
xlabel("Rw");
ylabel("Rq");
colormap winter;
colorbar;

subplot(2,3,3);
imagesc(Rw,Rq,F.q3);
%contourf(Rw,Rq,F.q3)
xlabel("Rw");
ylabel("Rq");
colormap winter;
colorbar;

subplot(2,3,4);
imagesc(Rw,Rq,F.wx);
%contourf(Rw,Rq,F.wx)
xlabel("Rw");
ylabel("Rq");
colormap winter;
colorbar;

subplot(2,3,5);
imagesc(Rw,Rq,F.wy);
%contourf(Rw,Rq,F.wy)
xlabel("Rw");
ylabel("Rq");
colormap winter;
colorbar;

subplot(2,3,6);
imagesc(Rw,Rq,F.wz);
%contourf(Rw,Rq,F.wz)
xlabel("Rw");
ylabel("Rq");
colormap winter;
colorbar;
