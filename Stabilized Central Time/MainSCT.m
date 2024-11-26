%% 
%% Question 4: Stabilized Centered Time


a = 0;
b = 1;
Tf = 0.5;
f = @(x) (x >= 0 & x <= 1) .* sin(pi*x) + (x < 0 | x > 1) .* 0;
u_Exact = @(x,t) exp(-pi^2*t) * sin(pi*x);


% % Error Plot
% 
N = [80,90,100];
dts = zeros(1,length(N));
% 
err_2 = zeros(1,length(N));
err_Inf = zeros(1,length(N));
% 
for i = 1:length(N)
    [dt,err2,err_inf] = SCT(a,b,N(i),Tf,f,u_Exact,i);
    err_2(i) = err2;
    err_Inf(i) = err_inf;
    dts(i) = dt;
end

% 
% Order of convergence

p2 = polyfit(log(dts),log(err_2),1);
disp(floor(p2(1)))
pInf = polyfit(log(dts),log(err_Inf),1);
disp(floor(pInf(1)))

figure(2*i+1)
loglog(dts, err_2, '-r*'); hold on;
loglog(dts, err_Inf, '-k*'); hold on;
xlabel('Final Time T_f');
ylabel('Error')
legend('Norm_2','Norm_{inf}')
title(sprintf('Convergence plot (Order of convergence = %.0f)', floor(p2(1))));
hold off;

%%




