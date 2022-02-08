theta = [2.82, 2.58, 2.59, 2.58, 2.87, 2.59, 2.59, 2.55, 2.83, 2.54, 2.48, 3.05, 2.74, 2.43, 2.71, 3.35];
N=[97, 98, 98, 98, 97, 98, 98, 98, 97, 98, 98, 96, 97, 98, 97, 95];

theta_m = [4.55, 3.78, 4.51, 2.18, 2.21, 2.15, 2.84, 2.56, 3.10, 2.36, 2.20, 3.29, 2.21, 2.78, 2.16, 2.54];
N_m = [89, 91, 89, 98, 98, 98, 95, 96, 94, 97, 98, 93, 98, 95, 98, 96, 97];

a = 16.63042;
b = 99.75063;

y = 0.0:0.01:18;

year = string([1893:1908])

v = zeros(length(theta), length(y));



% the areas under the curve (y,v).
areas_f = theta*0;

%legend(year)

subplot(1,2,1);

hold on;
for i = 1 : length(theta)
    for j = 1 : length(y)
      obj = @(x)func(x, theta(i), y(j), a, b);
      v(i,j) = integral(obj, 0.1, theta(i)/y(j), 'AbsTol', 1e-15, 'RelTol', 1e-5);
    end
    plot(y, v(i,:));
    tail_y = y;
    tail_y(1:((100-N(i))*100))=[];
    tail_v=v;
    tail_v(:,(1:((100-N(i))*100)))=[];
    areas_f(i) = trapz(tail_y, tail_v(i,:));
end

hold off;
title('Female')

v = zeros(length(theta), length(y));

subplot(1,2,2)

% the areas under the curve (y,v).
areas_m = theta*0;
hold on;
for i = 1 : length(theta_m)
    for j = 1 : length(y)
      obj = @(x)func(x, theta_m(i), y(j), a, b);
      v(i,j) = integral(obj, 0, theta_m(i)/y(j), 'AbsTol', 1e-15, 'RelTol', 1e-5);
    end
    plot(y, v(i,:));
    tail_y = y;
    tail_y(1:((100-N_m(i))*100))=[];
    tail_v=v;
    tail_v(:,(1:((100-N_m(i))*100)))=[];
    areas_m(i) = trapz(tail_y, tail_v(i,:));
end
hold off;
title('Male')

legend(year)


% t is "theta"
function r = func(x, t, y, a, b)
r = b^a/gamma(a)/t*exp(-b.*(x)).*(1-y/t*x).^(-1+1./x).*(x).^(a-1);
end