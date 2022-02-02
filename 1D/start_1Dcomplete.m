close all
clear
clc
N=50% number of shells 
k=round(0.1*N)+1 %k=round(0.1*N)+1; % k is off number of shells between NO^** and NO^+. set to approximately 10% of total shells 
n0=50;
dp=0.1;  %   peak density in um^-3 or 1e12 cm^-3
gamma1=100;
gamma2=100;

t_final=1e5;  % times are in ns 
dt=1e3;
tspan=[0:dt:t_final];

% r(theta)= a*b/sqrt(b^2 cos)

%     sigma_x=1*1000;
%     sigma_y=0.55*1000;  %
%     sigma_z=0.70*1000;  % Gaussian width convert from mm to um defualt 0.42


a=1000; % convert mm to nm
b=0.5*1000;

env=5;  % sigma environment

M=50;
d_phi=2*pi/M;

    sigma=b;
    phi=0;
    [t,rs, edens, n_ions, n_highs, nDens, Ts, Tn, r0, tau, ns]=complete1D(N,k,dp,n0,sigma,env,tspan,gamma1,gamma2, phi, d_phi);

 %%
 t=t/1000;

figure()
plot(t,rs)
ylabel('r')

figure()   
plot(t,n_ions)
ylabel('NO^+ ')
xlabel('us')

figure()
plot(t,n_highs)
ylabel('NO^{**} ')
xlabel('us')

figure()
plot(t,nDens)
ylabel('NO^*')
xlabel('us')

    