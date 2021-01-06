%[x,y] = ode4(@odefun,h, T, f0)
%odefun дифференциальное уравнение вида y' = f(x,y)
%h - шаг по сетке
%T = [x0, xf] x0 - нижний предел интегрирования, xf - верхний
%f0 = [x0,y0] - начальные условия
function [x,y] = ode4(odefun,h, T, f0)
x = T(1):h:T(2);
y = zeros(length(f0(2:end)),length(x));
y(:,1) = f0(2:end);
x(1) = f0(1);
    for n = 1:length(x)-1
        %%вычисление коэффициентов
        k1 = odefun(x(n), y(:,n));
        k2 = odefun(x(n)+h/2, y(:,n)+h*k1/2);
        k3 = odefun(x(n)+h/2, y(:,n)+h*k2/2);
        k4 = odefun(x(n)+h, y(:,n)+h*k3);
        %%вычисление приближенного значения по итерационной формуле
        y(:,n+1) = y(:,n)+h/6*(k1+2*k2+2*k3+k4);
    end
end