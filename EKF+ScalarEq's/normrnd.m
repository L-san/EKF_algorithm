function y = normrnd(a,sigma)
    sz = length(a);
    u = rand(12,sz);
    x=sum(u)'-6;
    y = sigma.*x+a;
end