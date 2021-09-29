function Q = quat_mult(a,b)
Q(1) = a(1)*b(1)-dot(a(2:4),b(2:4));
Q(2:4) = a(1)*b(2:4)+cross(a(2:4),b(2:4))+b(1)*a(2:4);
Q = Q';
end