clear variables; clc;
n = 10;

Rcq = 1e15:(1e18-1e15)/(n-1):1e18;
Rcw = 1e15:(1e18-1e15)/(n-1):1e18;

Qcq = 0:1e10/(n-1):1e10;
Qcw = 0:1e10/(n-1):1e10;

Pcq = 1e6:(1e15-1e6)/(n-1):1e15;
Pcw = 1e6:(1e15-1e6)/(n-1):1e15;

F.q1 = ones(n,n,n,n,n,n);
F.q2 = ones(n,n,n,n,n,n);
F.q3 = ones(n,n,n,n,n,n);
F.wx = ones(n,n,n,n,n,n);
F.wy = ones(n,n,n,n,n,n);
F.wz = ones(n,n,n,n,n,n);

for a = 1:n
    for b = 1:n
        R = diag([Rcq(a) Rcq(a) Rcq(a) Rcw(b) Rcw(b) Rcw(b)]);
        for c = 1:n
            for d = 1:n
                P = diag([Pcq(c) Pcq(c) Pcq(c) Pcw(d) Pcw(d) Pcw(d)]);
                for j = 1:n
                    for g = 1:n
                        Q = diag([Qcq(j) Qcq(j) Qcq(j) Qcw(g) Qcw(g) Qcw(g)]);
                        f = main(P,R,Q);
                        if(isnan(f))
                            F.q1(a,b,c,d,j,g) = 0;
                            F.q2(a,b,c,d,j,g) = 0;
                            F.q3(a,b,c,d,j,g) = 0;
                            F.wx(a,b,c,d,j,g) = 0;
                            F.wy(a,b,c,d,j,g) = 0;
                            F.wz(a,b,c,d,j,g) = 0;
                        else
                            F.q1(a,b,c,d,j,g) = f(1);
                            F.q2(a,b,c,d,j,g) = f(2);
                            F.q3(a,b,c,d,j,g) = f(3);
                            F.wx(a,b,c,d,j,g) = f(4);
                            F.wy(a,b,c,d,j,g) = f(5);
                            F.wz(a,b,c,d,j,g) = f(6);
                        end
                    end
                end
            end
        end
    end
end

