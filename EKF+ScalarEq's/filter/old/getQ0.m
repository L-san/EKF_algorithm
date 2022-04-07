function ans = getQ0(att)
ans = real(sqrt(1-norm(att)^2));
end