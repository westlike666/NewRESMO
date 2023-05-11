
function [t,rs, edens, n_ions, n_highs, nDens,Ts, Tn, r0, tau, ns]=complete1D(N,k,dp,n0,purpose,all_pqn,sigma,env,tspan,gamma1,gamma2, phi, d_phi)

%% define initial constant
kB = 1.3806504e-23;                 % #m2 kg s-2 K-1;J K-1; 8.62e-5 #eV K-1,
Ry = 2.179872e-18;                  % #J: Rydberg energy in J 
m_NO=30.006/1000/6.02214076e23;    % kg

kBm=kB/m_NO*1e12/1e18; %k_b / m_NO #um2 ns-2 K-1

%kBm=0.0002770924685030197; %k_b / m_NO #um2 ns-2 K-1

a0=5.2917721092e-5;     % bohr radius in um
%Rn0=n.^2*a0;            % radius of Rydb. atom by bohr model using semi-classical method

sigma_z=sigma;  % Gaussian width convert from mm to um.

T_e=50; T_i=0;   % K

tau=(1/kBm*sigma^2/(T_e+T_i))^0.5;  % ns  let sigma_0 equals to 1000um 

s0=1;
%tau=1;
%k=round(0.1*N)+1; % k is off number of shells between NO^** and NO^+. set to approximately 20% of total shells 
t1=(2^(2/k)-1)^0.5*tau; % t1/tau is always constant no matter of tau
s=zeros(N+1,1)';


for n=[1:length(s)]
    s(n)=s0.*(1+(t1/tau).^2).^((n-1)/2);  % defining the boudaries of initial shells according to s(t1)=s(0)*(1+t1^2/tau^2)^0.5, such that they all overlap at the same time
end

r_z=s/max(s); % normalize

n=[1:length(s)-1];
Tn=tau*((1+(t1/tau)^2).^(n)-1).^0.5;  % the overlapping times can be calculated by s_1(t_n)=s_n(t_1). ignoring Tn=0. There are n boudaries and n-1 Tn



%n0=50;
%env=5;  %sigma enviroment
r_z=r_z*env*sigma; % N+1 columes
%penning_fraction=0.45;

%% fill shells for NO^+, NO^*, NO^**
%dp=1; % peak density
%n0=50

den0=dp*exp(-([0,r_z(2:end-1)].^2)./(2.*sigma_z^2));  % 1*N row
%den0=1000*dp*exp(-([0,r_z(2:end-1)].^2)./(2.*sigma_z^2))./(2*pi)^0.5./sigma_z;  % 1*N row

N0=n0*ones(1,N);
[PF,eden,n_Ryd]=penningfraction(N0,den0); % density in um^-3
n_ion=eden;

% n_ion=penning_fraction/2*den0; % NO^+ ion 
% eden=n_ion; % electron 
% n_Ryd=(1-penning_fraction)*den0;  % remaining Rydberg that not involved in ionization 


n_high=0*n_ion; % long-lived NO^**

expand=all_pqn;  % expand the NO^* Rydberg into different levels 
if expand
     n_min=10;
     n_max=80;
     f=@(x)5.95*x.^5;                % This is the penning fraction distribution
     np=n_min:fix(n0/sqrt(2));      % Array of n states allowed after Penn ion
     ind=1:length(np); 
     nl=[n_min:n_max];
     ns=length(nl);
     n_dist=0*nl;
     n_dist(ind)=f(np/n0)/sum(f(np/n0));% dividing by the sum normalizes the function
     nDen=n_dist'*eden;
     nDen(nl==n0,:)=n_Ryd;
     %T_penning=(-Ry*den0./n0^2 + Ry*n_Ryd./n0^2 + Ry*sum(nDen(ind,:)./(nl(ind).^2)',1))./(3/2*kB.*eden); % by energy conservation, initial before penning is zero 
     T_penning=sum(-Ry*den0./n0^2 + Ry*n_Ryd./n0^2 + Ry*sum(nDen(ind,:)./(nl(ind).^2)',1))./sum(3/2*kB.*eden); % by energy conservation, initial before penning is zero 

     T_e=T_penning % initial temperature. comment if wish to use preset Temperature.
else
     nl=n0;
     ns=1;
     nDen=den0-n_ion;
end 

r0=r_z'; 

u0=0*r0;

y0=[r0; u0; eden'; n_ion'; n_high'; reshape(nDen,[ns*N,1]);T_e];

[ni,nf,II,minn,maxn,diffsn]=buildns(nl); 

% phi=pi/4;
% d_phi=pi/100;

%% rate equations
function dy=ode1D(t,y)

r=y(1:N+1);  m=N+1;% there are N+1 boudaries
U=y(m+1:m+N+1); m=m+N+1; % radial velocity
eden=y(m+1:m+N);  m=m+N; % Nx1 colume; electrons 
n_ion=y(m+1:m+N); m=m+N; % Nx1 colume;  NO^+ ion    
n_high=y(m+1:m+N);  m=m+N; % Nx1 colume;  long-lived NO^** Rydberg 
nDen=y(m+1:m+ns*N); m=m+ns*N; % Nx1 colume;  short-lived NO^* Rydberg
T=y(m+1); m=m+1; % scaler: electron temperature 


nDen=reshape(nDen,[ns,N]); % ns x N martrix or row

U=t/(t^2+tau^2).*r; % expanding volecity at r: N+1 rows
u=U(1:end-1);
rr=r(1:end-1);


dU=mean(-kBm*T*diff(eden)./eden(1:end-1)./diff(rr)./rr(1:end-1))*r;  % averaging gradient 

% du=-kBm*T*diff(eden)./eden(1:end-1)./diff(rr)./rr(1:end-1).*rr(1:end-1); 
% dU=[du(1);du; du(end)];
dU=0*dU;

%infintesmall volume dV=sin(phi) r^2 dr d_phi d_theta, integrate over dr
%from 0 to r and d_theta from 0 to 2*pi to find total volume at r(theta),
%the diff calculate the volume of slices betwenn shell boundaries


V=2*pi/3*d_phi*(diff(r.^3)); % N*1 colume volume of each shell, expand volume to 3D


if purpose=='soc' 
    dV=2*pi*d_phi*(diff(r.^2.*U));

elseif purpose=='geo' 
%    V=diff(r);
%    dV=diff(U);
    dV=0;  % This is to aviod the artefact of bifurcation in images 
end




k_ion=kION(nl,T); % 1* ns row.   ionization
k_tbr=kTBR(nl,T); % 1* ns row.   three body recombination

%k_n_np=zeros(ns,ns);
k_n_np=knnp(ni,nf,II,minn,maxn,diffsn,T);  % ns * ns matrix

k_pd=kPD(nl); % 1* ns row.  predissociation
k_DR=kDR(T); % scaler.  dissociative recombination 



% k_CT=gamma1*k_tbr.^(1/2).*u;  %  N*ns matrix. charge transfer rate: larger speed (r) causes more charge transfer within that shell at given amount of time;
% k_CT2=gamma2*max(k_tbr).^(1/2).*u; % N*1 colume, since NO^** is not distinguished by pqn 

% k_CT=gamma1.*diff(r).*(a0*nl.^2)/T;  %  N*ns matrix. charge transfer rate: larger speed (r) causes more charge transfer within that shell at given amount of time;
% k_CT2=gamma2.*diff(r).*(a0.^2); % N*1 colume, since NO^** is not distinguished by pqn 

k_CT=1e5*gamma1.*u.*(a0*nl.^2)/T;  %  N*ns matrix. charge transfer rate: larger speed (r) causes more charge transfer within that shell at given amount of time;
k_CT2=1e4*gamma2.*u.*(a0.^2); % N*1 colume, since NO^** is not distinguished by pqn 

 k_CT=t/(t^2+tau^2)*1e-5*gamma1.*rr;
 k_CT2=1e-5*gamma2.*rr;
%%



%e + NO^+ +NO^* == 2NO^**   with k_CT
%e + NO^+ +NO^** == 2NO^**  with k_CT2   

%e + NO^* ==NO^+ + 2e       k_ion
%2e + NO^+ ==NO^* + e       k_tbr


% Now consider the reactions happening in the overlapping shells of NO^+ (moving with u), NO^* (stationary) , NO^** (moving with 0.5u and falling k shells behind NO^+ )     
%This means that the firth kth NO^** shells are alaways empty. and the kth frontier NO^+ should not contribute to the NO^** production since there is no overlap with NO^** frontiers    
% The overlapping times are determined by Tn

step=sum(t>=Tn)+1;

d_ionize=k_ion.*eden.*[nDen(:,step:end),zeros(ns,step-1)]'; % N*ns matrix: the outer (step) shells are zero
% d_ionize(1:end-step+1,:) is the overlapping (non-zero) region

d_tbr=k_tbr.*eden.^4; % N*1 colume. one eden power is aborbed in k_tbr 
% d_tbr(1:end-step+1,:) is the overlapping (non-zero) region

d_n_np=sum(k_n_np,2).*nDen.*[zeros(1,step-1),eden(1:end-step+1)']; % ns*N matrix

%d_n_np=sum(k_n_np,2).*nden.*eden'; % ns*N matrix

k_np_n=k_n_np';

d_np_n=k_np_n*nDen.*[zeros(1,step-1),eden(1:end-step+1)'];      % ns*N matrix



d_CT=k_CT.*n_ion.*[nDen(:,step:end),zeros(ns,step-1)]'; % N*ns matrix: the outer (step) shells are zero
% d_CT(1:end-step+1,:) is the overlapping (non-zero) region

d_eden=sum(-d_CT+d_ionize-d_tbr,2); % N*1 colums 

d_ion=d_eden;

d_nDen=[zeros(step-1,ns);-d_CT(1:end-step+1,:)-d_ionize(1:end-step+1,:)+d_tbr(1:end-step+1,:)];  % N*ns matrix: the inner (step) shells are zero 

gate=1;

if step<k+1 % let the NO^** wait for k shells before starting move together with NO^+ shells. Within this time the NO^** align with NO^*
    d_high1=gate*[zeros(step-1,ns);2*d_CT(1:end-step+1,:)];  % N*ns matrix: the inner (step) shells are zero; or set all to zero   
    
    d_CT2= k_CT2.*n_ion.*[n_high(step:end);zeros(step-1,1)]; % N*1 colume
    
    % d_CT2(1:end-step+1) is the overlapping (non-zero) region
    
    
    d_high2=gate*[zeros(step-1,1);d_CT2(1:end-step+1)];  % N*1 colume
    
else
    d_high1=[zeros(k,ns);2*d_CT(1:end-step+1,:);zeros(step-1-k,ns)];  % N*ns matrix; NO^** start moving, falling behind of NO^+ by k shells. so the overall velocity is u/2.  
    
    d_CT2= k_CT2.*n_ion.*[n_high(k+1:end);zeros(k,1)]; % N*1 colume
    
    % d_CT2(1:end-k) is the overlapping (non-zero) region
    
    d_high2=[zeros(k,1);d_CT2(1:end-k)];  % N*1 colume   
end
 
% d_ion2=-d_CT2;
 d_eden2=-d_CT2;
 
 %D_eden_old=d_eden+d_eden2; % to calculate temperature before energy loss
 D_eden_old=sum(d_ionize-d_tbr,2);
 
 D_eden=d_eden+d_eden2-k_DR*n_ion.*eden-eden.*dV./V; % consider also the expansion of volume
 
 %d_eden=d_eden+d_eden2-k_DR*n_ion.*eden-eden.*dV./V; 
 D_ion=D_eden;
 
% D_nDen_old=d_nDen-d_n_np'+d_np_n'; % to calculate temperature before energy loss
 D_nDen_old=[zeros(step-1,ns);-d_ionize(1:end-step+1,:)+d_tbr(1:end-step+1,:)]-d_n_np'+d_np_n';  % N*ns matrix: the inner (step) shells are zero 
 
 %D_nDen=d_nDen-d_n_np'+d_np_n'-k_pd.*nDen'-nDen'.*dV./V; % assume Rydberg shells expand, 
 D_nDen=d_nDen-d_n_np'+d_np_n'-k_pd.*nDen'; % assume Rydberg not expand, 
 
 d_high=sum(d_high1,2)+d_high2; % N*1 colums, 

 D_high=d_high-n_high.*dV./V; % the 1/8 factor cancels in dV/V

 dT1=(Ry*sum(sum(D_nDen_old./nl.^2.*V))-1.5*kB*T*sum(D_eden_old.*V))./(1.5*kB*sum(eden.*V)); % the total energy conserve before the dissipation. Ry and kB both in S.I. unit 

 dT2=-2*t/(t^2+tau^2)*T; % temperature change due to expansion
 
%  if isreal(dT1)==0
%      dT1=0;
%  end
 
 dT=dT1+dT2;


dy=[U; dU; D_eden; D_ion; D_high; reshape(D_nDen',[N*ns,1]);dT];

end


%% ouputting the ode results


% t_final=max(Tn);
% 
 [t,y]=ode45(@(t,y)ode1D(t,y),tspan,y0);
% 

rs=y(:,1:N+1);  m=N+1; % N+1 boundaries of the shell 
Us=y(:,m+1:m+N+1); m=m+N+1; %Velocities 
edens=y(:,m+1:m+N);  m=m+N; % Nx1 colume; electrons 
n_ions=y(:,m+1:m+N); m=m+N; % Nx1 colume;  NO^+ ion    
n_highs=y(:,m+1:m+N);  m=m+N; % Nx1 colume;  long-lived NO^** Rydberg 
nDens=y(:,m+1:m+ns*N); m=m+ns*N; % Nx1 colume;  short-lived NO^* Rydberg
Ts=y(:,m+1); m=m+1;
%rs=rs(:,1:end-1);
%nDens=reshape(nDens,[ns,N]); % ns x N martrix or row






%% visualize

%visualize=0;
% if visualize
% y_ana=@(t,y0,tau) y0.*(1+(t/tau).^2).^0.5;
% 
% figure()
% plot(r0,'.')
% xlabel('shells')
% ylabel('initial r')
% 
% figure()
% plot(t,y_ana(t,r_z,tau))
% hold on
% %plot(t,rs,'bx')
% xlabel('t')
% ylabel('radius(t)')
% 
% 
% 
% figure()
% plot(t,n_ions)
% xlabel('t')
% ylabel('NO^+ ')
% 
% figure()
% plot(t,nDens)
% xlabel('t')
% ylabel('NO^* ')
% 
% figure()
% plot(t,n_highs)
% xlabel('t')
% ylabel('NO^{**} ')


% if visualize
% y_ana=@(t,y0,tau) y0.*(1+(t/tau).^2).^0.5;
% figure()
% hold on 
% dt=0.1*tau;
% t_final=10*tau;
% for t=[0:dt:t_final]
%     s=y_ana(t,r_z,tau);
%     clf
%    % ylim([0 5])
%     t_n=Tn(abs(t-Tn)<dt);   
%     for i=[1:length(s)]
%         yline(s(i))
%         yline(0.5*(s(i)),'red') % start with k shells off 
% %        yline(0.5*(Y(i)+Y0(i)),'red') % start together 
%         yline(r_z(i),'green')
%     end
%     if isempty(t_n)==0
%         title(['t=',num2str(t_n)])
%     end    
%     pause(0.1)
%    
% end
% hold off
% end

end
