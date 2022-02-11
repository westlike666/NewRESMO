% This code can produce the 2D results of the shell model:
% After setting up the parameters, hit Run and wait  

close all
clear
clc

%% Adjustable parmeters to play with

% First choose the purpose of running.  'soc' or 'geo'
% There is a trade-off of computation cost between searching SOC conditions
% and showing 2D bifurcation images 

%purpose='soc';
purpose='geo'

if purpose=='soc'
    all_pqn=true;  M=1;
    % all_pqn decides wether to consder all principal quanutum number within [30,80] 
    % M is the number of azimuth angles in the ellipse spreading from 0 to 2*pi: M=1 produce a 1D result
    % This configs means it only run the 1D model once (M=1), but will populate all the pqn from 30 to 80 during the calculation 
    % Use this set-up to quickly search the conditions of self-organized-critical without genrating bifurcation images.  
elseif purpose=='geo'
     all_pqn=false; M=50
    % This configs will run the model 50 times at 50 azimuth angles to construct 2D images
    % Considering all pricipal quantum number in [30,80] is not recommnded due to computation time  
    % Use this set-up to genrate bifurcation images.  
end 

N=100 % number of shells 


n0=50; % initial principal quanutm number  
dp=0.1;  %   peak density in um^-3 or 1e12 cm^-3
gamma1=0.01; % This is the magnitude of charge transfer rate: k_CT=gamma1.*dr.*(a0*n.^2); The dimesion is arbitary

save=false % whether so save images as movie




%% parameters that are suggested to be fixed 
gamma2=0; % For now ignoring the charge transfer with high Rydberg: NO^+ e + NO** -> 2NO**

cmax=0.012; % color axis limit 
t_final=1e5;  % evolution time in ns 
dt=1e3; % steps of time to demonstrate
tspan=[0:dt:t_final]; 

k=round(0.1*N)+1; %k=round(0.1*N)+1; % k is off number of shells between NO^** and NO^+. set to approximately 10% of total shells 

% r(theta)= a*b/sqrt(b^2 cos)

%     sigma_x=1*1000;
%     sigma_y=0.55*1000;  %
%     sigma_z=0.70*1000;  % Gaussian width convert from mm to um defualt 0.42


a=1000; % width of ellipse, convert mm to nm
b=0.5*1000; % height of ellipse. 

env=5;  % sigma environment of density distribution 


%filename=['2D/'  'den=' num2str(dp), ' n0=' num2str(n0) ' a=', num2str(a), ' b=', num2str(b),  ' gamma1=' num2str(gamma1), ' gamma2=' num2str(gamma2) ' k_pd=' num2str(k_pd)];


tot_angle=2*pi; % the ploar angle run from 0 to pi/2*pi, meaning the results representing the upper hemisphere/ shpere
Phi=linspace(0,tot_angle,M); 
d_phi=tot_angle/M; 
R=a.*b./sqrt(b.^2.*cos(Phi).^2+a.^2.*sin(Phi).^2);

X=zeros(length(tspan),N,M);
Y=zeros(length(tspan),N,M);
Rs=zeros(length(tspan),N+1,M);
N_ions=zeros(length(tspan),N,M);
N_highs=zeros(length(tspan),N,M);
NDens=zeros(length(tspan),N,M);



for m=[1:M]
    sigma=R(m);
    phi=Phi(m);
    [t,rs, edens, n_ions, n_highs, nDens, Ts, Tn, r0, tau, ns]=complete1D(N,k,dp,n0,all_pqn,sigma,env,tspan,gamma1,gamma2, phi, d_phi);
    [x,y] = pol2cart(Phi(m),rs(:,1:end-1));
    X(:,:,m)=x;
    Y(:,:,m)=y;
    Rs(:,:,m)=rs;
    N_ions(:,:,m)=n_ions;
    N_highs(:,:,m)=n_highs;    
    NDens(:,:,m)=sum(reshape(nDens,[],ns,N),2);     
    
%     N_ions(:,:,m)=n_ions./diff(rs,1,2);
%     N_highs(:,:,m)=n_highs./diff(rs,1,2);    
%     NDens(:,:,m)=sum(reshape(nDens./diff(rs,1,2),[],ns,N),2); 
    
end 
%%
 t=t/1000; % converting time from ns to um

 if purpose=='soc'
     Vs=2*pi/3*d_phi*(diff(Rs.^3,1,2)); 
 elseif purpose=='geo'
     Vs=1; 
 end

 tot_ions=sum(sum(N_ions.*Vs,2),3);
 tot_nDens=sum(sum(NDens.*Vs(1,:,:),2),3);
 tot_highs=sum(sum(N_highs.*Vs/8,2),3);  % 1/2 rs yields 1/8 Vs
 
figure()

subplot(2,2,1)
plot(t,tot_nDens)
ylabel('NO^*')
%xlabel('us')

subplot(2,2,2)
plot(t,tot_ions)
ylabel('NO^+ ')
%xlabel('us')

subplot(2,2,3)
plot(t,tot_highs)
ylabel('NO^{**} ')
xlabel('us')

figure()
subplot(2,1,1)
plot(t,rs)
ylabel('r/um')
xlabel('us')

subplot(2,1,2)
plot(t,Ts)
ylabel('T_e/K')
xlabel('us')

% if save
%      v=VideoWriter([filename '.avi']);
%      open(v);
% end

if purpose=='geo'
    figure()
    for i=[1:length(tspan)]
        T=tspan(i);
        xx=reshape(X(i,:,:),[N,M]);
        yy=reshape(Y(i,:,:),[N,M]);
        x0=reshape(X(1,:,:),[N,M]);
        y0=reshape(Y(1,:,:),[N,M]);
        
        nn_ions=reshape(N_ions(i,:,:),N,M);
        nn_highs=reshape(N_highs(i,:,:),N,M);
        nnDens=reshape(NDens(i,:,:),N,M);
        
    %     nn_ions=reshape(N_ions(i,:,:).*Vs(i,:,:),N,M);
    %     nn_highs=reshape(N_highs(i,:,:).*Vs(i,:,:)/8,N,M);
    %     nnDens=reshape(NDens(i,:,:).*Vs(1,:,:),N,M);
        
        subplot(2,2,1)
    %    colormap winter
        h=pcolor(x0,y0,nnDens);
        set(h, 'EdgeColor', 'none')
        cmap=colormap;
        set(gca,'Color',cmap(1,:))
     %   caxis([0 cmax])
        title('NO^* Rydberg')
        %axis(15000*[-1 1 -1 1])
        
        subplot(2,2,2)
    %    colormap turbo
        h=pcolor(xx,yy,nn_ions);
        set(h, 'EdgeColor', 'none')
        cmap=colormap;
        set(gca,'Color',cmap(1,:))
    %    caxis([0 cmax])
        axis(15000*[-1 1 -1 1])
        title('NO^+ ions')
        
        subplot(2,2,3)
    %    colormap parula
        h=pcolor(0.5*xx,0.5*yy,nn_highs);
        set(h, 'EdgeColor', 'none')
        cmap=colormap;
        set(gca,'Color',cmap(1,:))
    %    caxis([0 cmax])
        caxis
        ylabel('relative density')
        axis(10000*[-1 1 -1 1]) 
        title('NO^{**} long-lived Rydberg')
        xlabel('r/um')
        
        subplot(2,2,4)
    %    colormap parula
        h=pcolor(xx,yy,nn_ions+[nn_highs(k+1:N,:); zeros(k,M)]);   % plot the overlapping region 
        set(h, 'EdgeColor', 'none')
        cmap=colormap;
        set(gca,'Color',cmap(1,:))
        %caxis([0 cmax])
        ylabel('relative density')
        axis(15000*[-1 1 -1 1]) 
        title('NO^{**} and NO^+')        
        xlabel('r/um')
        
        if save
            frame=getframe(gcf);
            writeVideo(v,frame); 
        else
            pause(0.1)
        end
    end    

    if save
    close(v)
    end
end





