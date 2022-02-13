function Jh = jacobiH(attitude, Borb)
BAx = Borb(1); BAy = Borb(2); BAz = Borb(3);
q1 = attitude(1); q2 = attitude(2); q3 = attitude(3); q4 = attitude(4);

dbxdq1=2*(BAx*q1-BAz*q3+BAy*q4);
dbxdq2=2*(BAx*q2+BAy*q3+BAz*q4);
dbxdq3=-2*(BAz*q1-BAy*q2+BAx*q3);
dbxdq4=2*(BAy*q1+BAz*q2-BAx*q4);

dbydq1=2*(BAy*q1+BAz*q2-BAx*q4);
dbydq2=2*(BAz*q1-BAy*q2+BAx*q3);
dbydq3=2*(BAx*q2+BAy*q3+BAz*q4);
dbydq4=-2*(BAx*q1-BAz*q3+BAy*q4);

dbzdq1=2*(BAz*q1-BAy*q2+BAx*q4);
dbzdq2=-2*(BAy*q1+BAz*q2-BAx*q3);
dbzdq3=2*(BAx*q2-BAz*q3+BAy*q4);
dbzdq4=2*(BAx*q1+BAy*q3+BAz*q4);

dbdq = [dbxdq1 dbxdq2 dbxdq3 dbxdq4;
        dbydq1 dbydq2 dbydq3 dbydq4;
        dbzdq1 dbzdq2 dbzdq3 dbzdq4];

dbdw = zeros(3,3);

dwdq = zeros(3,4);

dwdw = eye(3,3);

Jh = [0 0 0 0 0 0 0;
      dbdq dbdw;
      dwdq dwdw];
end