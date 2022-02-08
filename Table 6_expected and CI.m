%Create confidence intervals by bootstraping theta, a, and b. 
%Assume a and b follow normal distributions.
% standard errors of theta: 
%Female:
%0.6296516 0.5468509 0.5550822 0.5704058 0.6897220 0.6061743 0.6148073 0.6211084 0.7698563 0.6519009 0.6623791 0.9968648 0.8349735 0.6901991 0.8523762 1.282152
%All female cohorts: 0.6240818
%Male:
%2.6690733 1.8622769 0.6399785 0.6382735 0.6279511 0.6396121 0.9935913 0.8599265 1.1412786 0.7347663 0.6365612 1.3093119 0.6364934 0.9929886 0.6408177 0.8682877
%All male cohorts: 
%0.7449179

SE_theta_F = [0.6296516 0.5468509 0.5550822 0.5704058 0.6897220 0.6061743 0.6148073 0.6211084 0.7698563 0.6519009 0.6623791 0.9968648 0.8349735 0.6901991 0.8523762 1.282152];
SE_theta_F_mean = 0.62;
SE_theta_M = [2.6690733 1.8622769 0.6399785 0.6382735 0.6279511 0.6396121 0.9935913 0.8599265 1.1412786 0.7347663 0.6365612 1.3093119 0.6364934 0.9929886 0.6408177 0.8682877];
SE_theta_M_mean = 0.74;
theta = [2.82, 2.58, 2.59, 2.58, 2.87, 2.59, 2.59, 2.55, 2.83, 2.54, 2.48, 3.05, 2.74, 2.43, 2.71, 3.35];
N=[97, 98, 98, 98, 97, 98, 98, 98, 97, 98, 98, 96, 97, 98, 97, 95];

theta_m = [4.55, 3.78, 4.51, 2.18, 2.21, 2.15, 2.84, 2.56, 3.10, 2.36, 2.20, 3.29, 2.21, 2.78, 2.16, 2.54];
N_m = [89, 91, 89, 98, 98, 98, 95, 96, 94, 97, 98, 93, 98, 95, 98, 96, 97];


%shape a, rate b. 
std_a = 10.84530; %var
std_b = 402.14671; %var
cov_ab = 65.05094; %cov

%a = 3.38;
%b = 51.75;

mu=[16.63042 99.75063];
Sigma =[std_a cov_ab; cov_ab std_b];
rng('default') %For reproducibility
rand_ab=mvnrnd(mu, Sigma, 50); %generated random a and b


age=[100, 105];

 % CI: the areas under the curve (y,v).
    areas_f = zeros(length(rand_ab), length(theta));
        % the areas under the curve (y,v).
    areas_m = zeros(length(rand_ab), length(theta_m));
    
     % CI: the expected numbers.
    CI_f = zeros(length(age), 2, length(theta));
    CI_m = zeros(length(age), 2, length(theta_m));
    
 % Point: the areas under the curve (y,v).
%     p_areas_f = zeros(length(age), length(theta));
%         % the areas under the curve (y,v).
%     p_areas_m = zeros(length(age), length(theta_m));
    
     % Point: expected numbers (y,v).
    %p_f = zeros(length(age), length(theta));
        % the areas under the curve (y,v).
    %p_m = zeros(length(age), length(theta_m));
    


F_N_pop = [1242, 923, 986, 1003, 1491, 1117, 1139, 1121, 1715, 1267, 1319, 2742, 1992, 1429, 2136, 4287];
M_N_pop = [4179, 2614, 4382, 309, 285, 284, 917, 650, 1262, 456, 262, 1951, 299, 1021, 314, 666];


        
        
 % Confidence Intervals

for l = 1:length(age)
    age(l)
    
%     a=mu(1);
%     b=mu(2);
% 
%         y = 0.0:0.01:25;
% 
%         year = string([1893:1908]);
% 
%         v = zeros(length(theta), length(y));
% 
% 
%         for i = 1 : length(theta)
%             for j = 1 : length(y)
%               obj = @(x)func(x, theta(i), y(j), a, b);
%               v(i,j) = integral(obj, 0, theta(i)/y(j), 'AbsTol', 1e-15, 'RelTol', 1e-5);
%             end
%             tail_y = y;
%             tail_y(1:((age(l)-N(i))*100))=[];
%             tail_v=v;
%             tail_v(:,(1:((age(l)-N(i))*100)))=[];
%             p_areas_f(l,i) = trapz(tail_y, tail_v(i,:));
%         end
% 
% 
%         v = zeros(length(theta), length(y));
% 
% 
%         for i = 1 : length(theta_m)
%             for j = 1 : length(y)
%               obj = @(x)func(x, theta_m(i), y(j), a, b);
%               v(i,j) = integral(obj, 0, theta_m(i)/y(j), 'AbsTol', 1e-15, 'RelTol', 1e-5);
%             end
%             tail_y = y;
%             tail_y(1:((age(l)-N_m(i))*100))=[];
%             tail_v=v;
%             tail_v(:,(1:((age(l)-N_m(i))*100)))=[];
%             p_areas_m(l,i) = trapz(tail_y, tail_v(i,:));
%         end
%         
%     p_f(l,:)=F_N_pop.*p_areas_f(l,:)
%     p_m(l,:)=M_N_pop.*p_areas_m(l,:)
    
        for k = 1: length(rand_ab)
            a=rand_ab(k, 1);
            b=rand_ab(k, 2);

            y = 0.0:0.01:25;

            year = string([1893:1908]);

            v = zeros(length(theta), length(y));

            for i = 1 : length(theta)
            theta_i = normrnd(theta(i), SE_theta_F_mean);
            while theta_i < 0.0001
              theta_i = normrnd(theta(i), SE_theta_F_mean);
            end
            for j = 1 : length(y)
              obj = @(x)func(x, theta_i, y(j), a, b);
              v(i,j) = integral(obj, 0, theta_i/y(j), 'AbsTol', 1e-12, 'RelTol', 1e-5);
            end
            tail_y = y;
            tail_y(1:((age(l)-N(i))*100))=[];
            tail_v=v;
            tail_v(:,(1:((age(l)-N(i))*100)))=[];
            areas_f(k, i) = trapz(tail_y, tail_v(i,:));
        end


            y = 0.0:0.01:25;

            year = string([1893:1908]);

            v = zeros(length(theta), length(y));
            for i = 1 : length(theta)
            theta_mi = normrnd(theta_m(i), SE_theta_M_mean);
            while theta_mi < 0.0001
              theta_mi = normrnd(theta_m(i), SE_theta_M_mean);
            end
            for j = 1 : length(y)
              obj = @(x)func(x, theta_mi, y(j), a, b);
              v(i,j) = integral(obj, 0, theta_mi/y(j), 'AbsTol', 1e-12, 'RelTol', 1e-5);
            end
            tail_y = y;
            tail_y(1:((age(l)-N_m(i))*100))=[];
            tail_v=v;
            tail_v(:,(1:((age(l)-N_m(i))*100)))=[];
            areas_m(k, i) = trapz(tail_y, tail_v(i,:));
        end

    end

    confidence_f = zeros(2, length(theta));
    confidence_f(1,:)=quantile(areas_f, 0.025);
    confidence_f(2,:)=quantile(areas_f, 0.975);

    confidence_m = zeros(2, length(theta));
    confidence_m(1,:)=quantile(areas_m, 0.025);
    confidence_m(2,:)=quantile(areas_m, 0.975);

    CI_f(l, :, :)=F_N_pop.*confidence_f
    CI_m(l, :, :)=M_N_pop.*confidence_m
    
    

end

% writematrix(p_f,'p_f.txt');
% writematrix(p_m,'p_m.txt');
writematrix(CI_f,'CI_f_14Oct.txt');
writematrix(CI_m,'CI_m_14Oct.txt');


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hold on;
% for k = 1:length(rand_ab)
%     a=rand_ab(k, 1);
%     b=rand_ab(k, 2);
%     y = 0.0:0.01:25;
%     i = 5;
%     j = 6;
%     obj = @(x)func(x, theta(i), y(j), a, b);
%     x = 0:0.01:3; %theta(i)/y(j);
%     target = obj(x);
%     plot(x,target);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%


% t is "theta"
function r = func(x, t, y, a, b)
r = b^a/gamma(a)/t*exp(-b.*(x)).*(1-y/t*x).^(-1+1./x).*(x).^(a-1);
end



