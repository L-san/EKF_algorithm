function [errN, err] = errcalc(x,y) %y - true %x - calc
n = height(x);
m = length(x);

%err = zeros(n,m);
% % 
% % Fy = 2*acos(sqrt(1-sum(y(1:3,:).^2)));
% % Fx = 2*acos(sqrt(1-sum(x(1:3,:).^2)));
% 
% err(1,:) = abs(Fx-Fy);
% 
% 
% if(~isnan(x))
%     for i = 1:n
%         for j = 1:m
%             err(i,j) = abs(y(i,j)-x(i,j));
%         end
%     end
%     errN = median(err')*180/pi;
% else
%     errN = ones(n,1)*NaN;
% end
% end

err = zeros(n,m);
if(~isnan(x))
    for i = 1:n
        for j = 1:m
            err(i,j) = abs(y(i,j)-x(i,j));
        end
    end
    errN = median(err');
else
    errN = ones(n,1)*NaN;
end
end