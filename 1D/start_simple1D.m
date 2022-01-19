close all
clear
clc
N=10% number of shells 
dp=1 % peak density
sigma=1 %
gamma1=50
gamma2=50


%t_final=10
visualize=true

t_final=10
dt=0.1
tspan=[0:dt:t_final]

[t,rs, edens, n_ions, n_highs, nDens, Tn, r0, tau]=simple1D(N,dp,sigma,tspan,gamma1,gamma2,visualize);


%t_final=max(Tn)
% figure()
% hold on 
% dt=0.1

% for j=[1:length(t)]
% %   s=y_ana(T,r0,tau);
%     T=t(j);
%     s=rs(j,:);  
%     clf
%     ylim([0 10])
%     t_n=Tn(abs(T-Tn)<dt);   
%     for i=[1:length(s)]
%         yline(s(i))
%         yline(0.5*(s(i)),'red') % start with k shells off 
% %        yline(0.5*(Y(i)+Y0(i)),'red') % start together 
%         yline(r0(i),'green')
%     end
%     if isempty(t_n)==0
%         title(['t=',num2str(t_n)])
%     end    
%     pause(0.1)
%    
% end

%% make video
filename=['N=' num2str(N) ' den' num2str(dp), 'sigma=', num2str(sigma) 'gamma1=' num2str(gamma1), 'gamma2=' num2str(gamma2)];
%dt=0.1;

y_ana=@(t,y0,tau) y0.*(1+(t/tau).^2).^0.5;
v=VideoWriter([filename '.mov']);
open(v);

rs=rs/sigma;
r0=r0/sigma;

figure()
for T=tspan
    s=y_ana(T,r0,tau);
    clf
  
    subplot(2,2,1)
    ylim([0 2])
%   t_n=Tn(abs(T-Tn)<=dt);
    t_n=T;
    for i=[1:length(s)]
        yline(s(i),'b--', 'linewidth',1)
        yline(0.5*(s(i)),'red', 'linewidth',1) % start with k shells off 
%        yline(0.5*(Y(i)+Y0(i)),'red') % start together 
        yline(r0(i),'green','linewidth', 1)
    end
%     if isempty(t_n)==0
%         title(['overlap t=',num2str(t_n) ' \tau_{exp}'])
%     end
    
    subplot(2,2,2)
    s_T=y_ana([0:dt:T],r0,tau);
    plot([0:dt:T],s_T)
    xlim([0 t_final])
    ylim([0 max(max(rs))])
    xlabel('t/ \tau_{exp}')
    ylabel('r/ \sigma_0')
    %pause(0.1)

    
    subplot(2,2,[3,4])
    xlabel('r/ \sigma_0')    
    xlim([0,y_ana(max(Tn),max(r0),tau)])
    xlim([0,10])
    ylabel('1D density')
    ylim([0,10])
    
%    index=(abs(T-t)<=2*dt);
    index=(T==tspan);  
    r_t=mean(rs(index,:),1);
    dr_t=[0,diff(r_t)];
    dr0=[0,diff(r0')];
    hold on
    for i=[1:length(s)]        
        s_nDen=mean(nDens(index,i),1);   
        bar(r0(i)-dr0(i)/2, s_nDen/dr0(i), dr0(i),'facecolor', 'green','edgecolor','none', 'facealpha',0.4) % devide the width to normalize 
        s_ion=mean(n_ions(index,i),1);   
        bar(r_t(i)-dr_t(i)/2, s_ion/dr_t(i), dr_t(i),'facecolor', 'blue','edgecolor','none', 'facealpha',0.4)
        s_high=mean(n_highs(index,i),1);   
        bar((r_t(i)-dr_t(i)/2)/2, s_high*2/dr_t(i), dr_t(i)/2,'facecolor', 'r','edgecolor','none', 'facealpha',0.8)
    end
    legend('NO^* Ryd',  'NO^+ ions', 'long-lived NO^{**}')
    hold off
    
    frame=getframe(gcf);
    writeVideo(v,frame); 
end

close(v)



