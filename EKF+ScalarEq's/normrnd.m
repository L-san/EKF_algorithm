function y = normrnd(a,sigma,size)
n = size(1);
m = size(2);
y = zeros(n,m);
for i = 1:m
    u = rand(12,n);
    x=sum(u)'-6;
    y(:,i) = sigma.*x+a;
end
end