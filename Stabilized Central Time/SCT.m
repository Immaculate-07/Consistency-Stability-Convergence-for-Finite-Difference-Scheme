function [dt,err2,err_inf] = SCT(a,b,N,Tf,f,u_Exact,i)
% Stabilaized Central Time.
% input: f is the initial condition function at t = 0 
%        a is the starting point of the spatial domain
%        b is the ending point of the spatial domain.
%        N: The number of time steps.
%        Tf: The final time.
%        u_Exact: is the exact solution
%        i is the index of the figure greater than 1
% output: dt is the time interval
%         err2 is the Norm 2 error
%         err_inf is the infinity norm
%
% Usage: [dt,err2,err_inf] = SCT(a,b,N(i),Tf,f,u_Exact,i) where f can be a compactly spp[0,1]


format long
% Step1: Discretize the meshes
M = 11; % Number of points in Spacial domain
dx = (b-a)/(M+1); % 
xi = linspace(a, b,M+2); % points in x axis
dt = Tf/N; 
tn = 0:dt:Tf; % points in t axis 


% Step2: Define the constant
lambda = (2*dt)/(dx^2);
alpha = lambda / (1+lambda);
beta = (1-lambda)/(1+lambda);


% Step3: Create the initial condition
U = zeros(M+2,length(tn)); % M - Spacial, N - Time


U(:,1) = f(xi); % Initial condition at t = 0;

UM1 = zeros(M,1); % Initialzize U^{-1}

% Step4: Set up the matrix A and B;
% Matrix A will be of dimension (M-2), since U(0,n)=U(M,n)=0 for all n

% Step I: Create the maindiagonal and offdiagonal elements

main_diag_A = 0 * ones(M,1) ; % Creating the diagonal vectors
off_diag_A = alpha * ones(M-1,1); % Creating the off diagonal vectors

% Step II: Create the matrix A and B;
 
A = diag(main_diag_A) + diag(off_diag_A,1) + diag(off_diag_A,-1); 


%Step5: Create the algorithm for each time interval;

for n = 1:length(tn)-1
    U(2:M+1, n+1) = A * U(2:M+1,n) + beta * UM1;
    UM1 = U(2:M+1,n);
    U(1,n+1) = 0; 
    U(M+2,n+1) = 0; 
end


% Display the Surface Plot
figure(2*i-1)
[X, T] = meshgrid(xi, tn);
surf(X, T, U');
xlabel('Space (x)');
ylabel('Time (t)');
zlabel('Temperature (U)');
title(sprintf('Foward Euler Solution of Heat Equation with dt=%.4f',dt));
    
%Step6 Get the solution(Exact and Numerical)
U_xi_Tf = U(:,end); % Approximated Solution at the final time

u_xi_Tf = arrayfun(@(x) u_Exact(x,Tf), xi); 

err2 = norm(U_xi_Tf'-u_xi_Tf,2);
err_inf = norm(U_xi_Tf'-u_xi_Tf,inf);

figure(2*i)
plot(xi,u_xi_Tf,'-r','LineWidth',2);
hold on;
plot(xi,U_xi_Tf,'b','LineWidth',2);
xlabel('x');
ylabel('U(x,t)');
title(sprintf('Plot of U_{Aprox} and u_{Exact} at T=%.2f',Tf));
legend('u_{Exact}', 'U_{Approx}');

hold off;
 
 
end

