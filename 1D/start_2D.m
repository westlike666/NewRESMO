close all
clear
clc
N=100% number of shells 
k=round(0.1*N)+1 %k=round(0.1*N)+1; % k is off number of shells between NO^** and NO^+. set to approximately 10% of total shells 
n0=50;
dp=0.5% peak density
gamma1=1
gamma2=10
k_pd=0.01;
k_DR=0.1

visualize=false;
save=false

cmax=0.012
t_final=700
dt=10
tspan=[0:dt:t_final];

% r(theta)= a*b/sqrt(b^2 cos)

%assume a>b
a=2;
b=1;

M=100;


filename=['2D/'  'den=' num2str(dp), ' n0=' num2str(n0) ' a=', num2str(a), ' b=', num2str(b),  ' gamma1=' num2str(gamma1), ' gamma2=' num2str(gamma2) ' k_pd=' num2str(k_pd)];


theta=linspace(0,2*pi,M);

R=a.*b./sqrt(b.^2.*cos(theta).^2+a.^2.*sin(theta).^2);

X=zeros(length(tspan),N,M);
Y=zeros(length(tspan),N,M);
N_ions=zeros(length(tspan),N,M);
N_highs=zeros(length(tspan),N,M);
NDens=zeros(length(tspan),N,M);



for m=[1:M]
    sigma=R(m);
    [t,rs, edens, n_ions, n_highs, nDens, Tn, r0, tau]=simple1D(N,k,dp,n0,sigma,tspan, gamma1,gamma2,k_pd,k_DR,visualize);
    [x,y] = pol2cart(theta(m),rs);
    X(:,:,m)=x;
    Y(:,:,m)=y;

    N_ions(:,:,m)=n_ions;
    N_highs(:,:,m)=n_highs;
    NDens(:,:,m)=nDens;
    
end 

if save
     v=VideoWriter([filename '.avi']);
     open(v);
end

for i=[1:length(tspan)]
    T=tspan(i);
    xx=reshape(X(i,:,:),[N,M]);
    yy=reshape(Y(i,:,:),[N,M]);
    x0=reshape(X(1,:,:),[N,M]);
    y0=reshape(Y(1,:,:),[N,M]);
    
    nn_ions=reshape(N_ions(i,:,:),N,M);
    nn_highs=reshape(N_highs(i,:,:),N,M);
    nnDens=reshape(NDens(i,:,:),N,M);
    
    subplot(2,2,1)
%    colormap winter
    h=pcolor(x0,y0,nnDens);
    set(h, 'EdgeColor', 'none')
    cmap=colormap;
    set(gca,'Color',cmap(1,:))
%    caxis([0 cmax])
    title('NO^* Rydberg')
    axis(10*[-1 1 -1 1])
    
    subplot(2,2,2)
%    colormap turbo
    h=pcolor(xx,yy,nn_ions);
    set(h, 'EdgeColor', 'none')
    cmap=colormap;
    set(gca,'Color',cmap(1,:))
%    caxis([0 cmax])
    axis(2500*[-1 1 -1 1])
    title('NO^+ ions')
    
    subplot(2,2,3)
%    colormap parula
    h=pcolor(0.5*xx,0.5*yy,nn_highs);
    set(h, 'EdgeColor', 'none')
    cmap=colormap;
    set(gca,'Color',cmap(1,:))
    caxis([0 cmax])
    ylabel('relative density')
    axis(2500*[-1 1 -1 1]) 
    title('NO^{**} long-lived Rydberg')
    xlabel('r/ \sigma_0')
    
    subplot(2,2,4)
%    colormap parula
    h=pcolor(xx,yy,nn_ions+[nn_highs(k+1:N,:); zeros(k,M)]);   % plot the overlapping region 
    set(h, 'EdgeColor', 'none')
    cmap=colormap;
    set(gca,'Color',cmap(1,:))
%    caxis([0 cmax])
    ylabel('relative density')
    axis(2500*[-1 1 -1 1]) 
    title('NO^{**} and NO^+')    
    
    xlabel('r/ \sigma_0')
    
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





