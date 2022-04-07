prec = 3*1e-3;
x = [1e5 1e5 0 0 0 0]';
n = 10000;
A = 1e2;
C1 = A/(1-exp(-n));
C2 = A-C1;
F = eye(6,1);
for i =1:n
    mu = C2+C1*exp(-i);
    P = diag([x(1) x(1) x(1) x(2) x(2) x(2)]);
    R = diag([x(3) x(3) x(3) x(4) x(4) x(4)]);
    Q = diag([x(5) x(5) x(5) x(6) x(6) x(6)]);
    f = main(P,R,Q);
    if(isnan(f))
        continue;
    else
        F = f';
    end
    Xt = x-mu*gradient(f)'/norm(gradient(f));
    mu = sqrt(n-1)-sqrt(i-1);
    if(abs(f)<prec)
        F = f'
        break;
    else
        x = Xt;
    end
end