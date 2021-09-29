function Q = quat_conj(a)
Q(1) = a(1);
Q(2:4) = -a(2:4);
Q = Q';
end