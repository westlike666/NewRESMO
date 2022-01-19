clc; clear; close all


%% constants 

kB = 1.3806504e-23;                 % #m2 kg s-2 K-1;J K-1; 8.62e-5 #eV K-1,
Ry = 2.179872e-18;                  % #J: Rydberg energy in J 
m_NO=30.006/1000/6.02214076e23;    % kg

kBm=kB/m_NO*1e12/1e18; %k_b / m_NO #um2 ns-2 K-1

kBm=0.0002770924685030197; %k_b / m_NO #um2 ns-2 K-1

sigma_z=0.7*1000;  % Gaussian width convert from mm to um.

T_e=5; T_i=0;   % K

tau=(1/kBm*sigma_z^2/(T_e+T_i))^0.5;  % ns

t=0

gamma=t/(tau^2+t^2);

%% define shells 
sigma_env=5;
N=10
r_z=linspace(sigma_env*sigma_z/N,sigma_env*sigma_z,N); %first define even posistions, not shells
s_z=exp([1:N]); % such that each shell experineces same amount of time of expansion: u(r,t)=gamma(t)* r 

s_z=s_z*sigma_env*sigma_z/max(s_z);% rescale shells such that they spread the whole sigma environment
d_s=diff(s_z);
plot(d_s./s_z(1:end-1)) % check the spacing of shell satisfies u(r,t)=gamma(t)* r; In another word, ds/s is a constant at time t

figure()
%plot(r_z,s_z)
%% defines quantities of NO^+, NO^* and NO^** in seperate shells

n0=80;
d_p=0.01;  % peak initial density at center 

den0=d_p.*exp(-(s_z.^2)./(2.*sigma_z^2)); %Gaussian distribution of density; devided into N shells

firstn=1;
numlev=100;                         % this is the number of n levels initially considered
deac=0;                             % start with all Rydberg states
nl=(firstn:firstn+numlev-1)';       % array of accessible n levels, from 1 to 100

ns=length(nl);                      %number of n-states

NL=zeros(ns,N);                     %matrix of shells per n level
N0=n0*ones(1,N);
[PF,eDen,rDen]=penningfraction(N0,den0);

 f=@(x)5.95*x.^5;                % This is the penning fraction distribution
 np=firstn:fix(n0/sqrt(2));      % Array of n states allowed after Penn ion
 nDen=nl*0;                      % This is the distribution of penning fraction of penning partner (Rydberg molecules) 
 ind=1:length(np); 
 nDen(ind)=f(np/n0)/sum(f(np/n0));% dividing by the sum normalizes the function
 NDEN=nDen*eDen;                  % eDen is equal to the number of the remaining penning partner
 NDEN(nl==n0,:)=rDen;             % set n0 to rden which is the remaining Rydberg molecules in state n0

 
 






