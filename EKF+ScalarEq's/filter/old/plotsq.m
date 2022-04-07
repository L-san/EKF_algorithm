
for i = 1:10
    figure;
    for j = 1:10
        subplot(2,5,j);
        imagesc(Rcq,Qcq,F.q1(1,1,:,:,1,1));
        contourf(Rcq,Qcq,F.q1(1,1,:,:,1,1))
        colormap winter;
        colorbar;
    end
end
%%
figure(1);
subplot(2,3,1);
surf(Rcq,Pcq,Qcq,F.q1(:,:,1,1,1,1));
xlabel("Rcq");
ylabel("Pcq");
zlabel("Qcq");
colormap winter;
colorbar;

subplot(2,3,2);
surf(Rcq,Pcq,Qcq,F.q2);
xlabel("Rcq");
ylabel("Pcq");
zlabel("Qcq");
colormap winter;
colorbar;

subplot(2,3,3);
surf(Rcq,Pcq,Qcq,F.q3);
xlabel("Rcq");
ylabel("Pcq");
zlabel("Qcq");
colormap winter;
colorbar;

subplot(2,3,4);
surf(Rcq,Pcq,Qcq,F.wx);
xlabel("Rcq");
ylabel("Pcq");
zlabel("Qcq");
colormap winter;
colorbar;

subplot(2,3,5);
surf(Rcq,Pcq,Qcq,F.wy);
xlabel("Rcq");
ylabel("Pcq");
zlabel("Qcq");
colormap winter;
colorbar;

subplot(2,3,6);
surf(Rcq,Pcq,Qcq,F.wz);
xlabel("Rcq");
ylabel("Pcq");
zlabel("Qcq");
colormap winter;
colorbar;
%%
%%
figure(2);
subplot(2,3,1);
surf(Rcw,Pcw,Qcw,F.q1);
xlabel("Rcw");
ylabel("Pcw");
zlabel("Qcw");
colormap winter;
colorbar;

subplot(2,3,2);
surf(Rcw,Pcw,Qcw,F.q2);
xlabel("Rcw");
ylabel("Pcw");
zlabel("Qcw");
colormap winter;
colorbar;

subplot(2,3,3);
surf(Rcw,Pcw,Qcw,F.q3);
xlabel("Rcw");
ylabel("Pcw");
zlabel("Qcw");
colormap winter;
colorbar;

subplot(2,3,4);
surf(Rcw,Pcw,Qcw,F.wx);
xlabel("Rcw");
ylabel("Pcw");
zlabel("Qcw");
colormap winter;
colorbar;

subplot(2,3,5);
surf(Rcw,Pcw,Qcw,F.wy);
xlabel("Rcw");
ylabel("Pcw");
zlabel("Qcw");
colormap winter;
colorbar;

subplot(2,3,6);
surf(Rcw,Pcw,Qcw,F.wz);
xlabel("Rcw");
ylabel("Pcw");
zlabel("Qcw");
colormap winter;
colorbar;

surf(Rcq,Pcq,F.q1(:,:,1,1,1,1));
contourf(Rcq,Pcq,F.q1(:,:,1,1,1,1))
xlabel("Rc");
ylabel("Pc");
colormap winter;
colorbar;