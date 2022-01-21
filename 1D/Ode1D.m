
close all
clear 
clc

%% define self-similar expanding shells  
N=15; % number of shells 

r0=1;
tau=1;
k=round(0.1*N)+1 % k is off number of shells between NO^** and NO^+. set to approximately 10% of total shells 
t1=(2^(2/k)-1)^0.5*tau;
%t1=1;
s=zeros(N,1)';


for n=[1:length(s)]
    s(n)=r0.*(1+(t1/tau).^2).^((n-1)/2);  % defining the initial shells according to s(t1)=s(0)*(1+t1^2/tau^2)^0.5, such that they all overlap at the same time
end

s0=s/max(s);

n=[1:length(s)];
Tn=tau*((1+(t1/tau)^2).^(n)-1).^0.5;  % the overlapping times can be calculated by s_1(t_n)=s_n(t_1). ignoring Tn=0.


%% fill shells for NO^+, NO^*, NO^**
n0=50;
sigma_z=1;
env=5;  %sigma enviroment
r_z=s0*env;
penning_fraction=0.45;
dp=1; % peak density
den0=dp*exp(-([0,r_z(2:end)].^2)./(2.*sigma_z^2)); 

N0=n0*ones(1,N);
[PF,eden,n_Ryd]=penningfraction(N0,den0);
n_ion=eden;



% n_ion=penning_fraction/2*den0; % NO^+ ion 
% eden=n_ion; % electron 
% n_Ryd=(1-penning_fraction)*den0;  % remaining Rydberg that not involved in ionization 

n_high=1*n_ion; % long-lived NO^**

expand=false  % expand the NO^* Rydberg into different levels 
if expand
     f=@(x)5.95*x.^5;                % This is the penning fraction distribution
     np=1:fix(n0/sqrt(2));      % Array of n states allowed after Penn ion
     ind=1:length(np); 
     nl=[1:100];
     ns=length(nl);
     n_dist=0*nl;
     n_dist(ind)=f(np/n0)/sum(f(np/n0));% dividing by the sum normalizes the function
     nDen=n_dist'*eden;
     nDen(nl==n0,:)=n_Ryd;
else
     ns=1;
     nDen=den0-n_ion;
end 


y0=[r_z'; eden'; n_ion'; n_high'; reshape(nDen,[ns*N,1])];


%% rate equations
t=6
y=y0    
r=y(1:N);  m=N; % let r represent the outer boundaries of the shell 
eden=y(m+1:m+N);  m=m+N; % Nx1 colume; electrons 
n_ion=y(m+1:m+N); m=m+N; % Nx1 colume;  NO^+ ion    
n_high=y(m+1:m+N);  m=m+N; % Nx1 colume;  long-lived NO^** Rydberg 
nDen=y(m+1:m+ns*N); m=m+ns*N; % Nx1 colume;  short-lived NO^* Rydberg

nDen=reshape(nDen,[ns,N]); % ns x N martrix or row

gamma_CT=1;

u=t/(t^2+tau^2).*r; % expanding volecity at r

k_CT=gamma_CT*u;  % calculate the charge transfer rate: larger speed (r) causes more charge transfer within that shell at given amount of time 
k_CT2=k_CT;

% Now consider the reactions happening in the overlapping shells of NO^+ (moving with u), NO^* (stationary) , NO^** (moving with 0.5u and falling k shells behind NO^+ )     
%This means that the firth kth NO^** shells are alaways empty. and the kth frontier NO^+ should not contribute to the NO^** production since there is no overlap with NO^** frontiers    
% The overlapping times are determined by Tn

step=sum(t>=Tn)+1

%e + NO^+ +NO^* == 2NO^**   with k_CT
%e + NO^+ +NO^** == 2NO^**  with k_CT2

d_CT=k_CT.*n_ion.*eden.*[nDen(:,step:end),zeros(ns,step-1)]'; % N*ns matrix: the outer (step) shells are zero

% d_CT(1:end-step+1,:) is the overlapping (non-zero) region

d_eden=sum(-d_CT,2); % N*1 colums 

d_ion=d_eden;

d_nDen=[zeros(step-1,ns);-d_CT(1:end-step+1,:)];  % N*ns matrix: the inner (step) shells are zero 

gate=1;

if step<k+1 % let the NO^** wait for k shells before starting move together with NO^+ shells. Within this time the NO^** align with NO^*
    d_high=gate*[zeros(step-1,ns);2*d_CT(1:end-step+1,:)];  % N*ns matrix: the inner (step) shells are zero; or set all to zero   
    
    d_CT2= k_CT2.*n_ion.*eden.*[n_high(step:end);zeros(step-1,1)]; % N*1 colume
    
    % d_CT2(1:end-step+1) is the overlapping (non-zero) region
    
    
    d_high2=gate*[zeros(step-1,1);d_CT2(1:end-step+1)];  % N*1 colume
    
else
    d_high=[zeros(k,ns);2*d_CT(1:end-step+1,:);zeros(step-1-k,ns)];  % N*ns matrix; NO^** start moving, falling behind of NO^+ by k shells. so the overall velocity is u/2.  
    
    d_CT2= k_CT2.*n_ion.*eden.*[n_high(1:end-k);zeros(k,1)]; % N*1 colume
    
    % d_CT2(1:end-k) is the overlapping (non-zero) region
    
    d_high2=[zeros(k,1);d_CT2(1:end-k)];  % N*1 colume   
end
 
 d_ion2=-d_CT2;
 d_eden2=-d_CT2;
 
 
 d_ion=d_ion+d_ion2;
 d_eden=d_eden+d_eden2;
 
d_high=sum(d_high,2)+d_high2; % N*1 colums


dy=[u; d_eden; d_ion; d_high; reshape(d_nDen',[N*ns,1])];




%% visualize
figure()
[t,y_int]=ode45(@(t,y)simpleode(t,y,tau),[0 10],s0);

plot(t,y_int)

figure()
y_ana=@(t,y0,tau) y0.*(1+(t/tau).^2).^0.5;

plot(t,y_ana(t,s0,tau))

figure()
plot(s0,'.')

figure()
hold on 
dt=0.1;
for t=[0:dt:10]
    s=y_ana(t,s0,tau);
    clf
    ylim([0 10])
    t_n=Tn(abs(t-Tn)<dt);   
    for i=[1:length(s)]
        yline(s(i))
        yline(0.5*(s(i)),'red') % start with k shells off 
%        yline(0.5*(s(i)+s0(i)),'red') % start together 
        yline(s0(i),'green')
    end
    if isempty(t_n)==0
        title(['t=',num2str(t_n)])
    end    
    pause(0.1)
   
end
hold off


function dy=simpleode(t,y,tau)
gamma=t/(t^2+tau^2);
dy=gamma.*y;
end


