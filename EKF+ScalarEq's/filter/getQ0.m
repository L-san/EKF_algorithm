function ans = getQ0(att)
ans = real(sqrt(1-sum(att.^2)));
end