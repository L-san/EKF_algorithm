function [errN, err] = errcalc(x,y)
n = height(x);
m = length(x);
err = zeros(n,m);
if(~isnan(x))
    for i = 1:n
        for j = 1:m
            err(i,j) = abs(1-x(i,j)/y(i+1,j))*100;
            if(isnan(err(i,j))||isinf(err(i,j)))
                err(i,j) = abs(x(i,j));
            end
            if(err(i,j)>100)
                err(i,j) = err(i,j-1);
            end
        end
    end
else
    errN = NaN;
end
errN = median(err');
end