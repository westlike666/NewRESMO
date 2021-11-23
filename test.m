%% initializing 

den0=[1 1];
n0=49;
vectorize=true

N=length(den0); %number of shells

% Set constants and lower cut-off for PQN distribution:
kB = 1.3806504e-23;                 % #m2 kg s-2 K-1;J K-1; 8.62e-5 #eV K-1,
Ry = 2.179872e-18;                  % #J: Rydberg energy in J 

kBm=0.0002770924685030197; %k_b / m_NO #um2 ns-2 K-1

firstn=1;
n_min=10;
numlev=100;                         % this is the number of n levels initially considered
deac=0;                             % start with all Rydberg states
nl=(firstn:firstn+numlev-1)';       % array of accessible n levels, from 1 to 100

ns=length(nl);                      %number of n-states

NL=zeros(ns,N);                     %matrix of shells per n level

DEN0=zeros(ns,N);
NDEN=NL*0;
EDEN=zeros(N,1);                    %electron densities at each shell
DEAC=zeros(N,1);                    
DEAC_DR=zeros(N,1);                 % N + O from NO^+ + e^- at each shell 
DEAC_PD=zeros(ns,N);                % N + O from NO^* of each accessible level at each shell 
DEAC_N_MIN=zeros(n_min,N);
T_PENNING=zeros(1,N);               %temperature at each shell after penning ionization


% for ii=1:N %loop though all shells
%     h=1+(ii-1)*ns;
%     k=ii*ns;
%     hh=1+(ii-1)*n_min;
%     kk=ii*n_min;
%     
%     
%     [PF,eden,rden]=penningfraction(n0,den0(ii));
%     % Redistributes the Penning partners over lower n's:
%     f=@(x)5.95*x.^5;                % This is the penning fraction distribution
%     np=firstn:fix(n0/sqrt(2));      % Array of n states allowed after Penn ion
%     ind=1:length(np);               % This is the distribution of penning fraction
%     nden=nl*0;
%     nden(ind)=eden*f(np/n0)/sum(f(np/n0));  % dividing by the sum normalizes the function
% 
%     nden(nl==n0)=rden;              % set n0 to rden
%     
%     % Set initial temperature:    (Robicheaux 2005 JPhysB)
%     T_PENNING(ii)=(-Ry*den0(ii)/n0^2 + Ry*rden/n0^2 + Ry*sum(nden(ind)./nl(ind).^2) )*1/(3/2*kB*eden); % by energy conservation
% 
%     deac=sum(nden(1:n_min));    % allow n<=n_min to decay
%     D_DEAC_N_MIN(:,ii)=nden(1:n_min); %
%     nden(1:n_min)=zeros(n_min,1); 
%     
%     NL(:,ii)=nl;       % save values for this shell in arrays
%     DEN0(:,ii)=repmat(den0(ii),ns,1);
%     NDEN(:,ii)=nden;     
%     EDEN(ii)=eden;
%     DEAC(ii)=deac;
% end
% 
% 
% DEN0
% NDEN
% EDEN
% DEAC
% D_DEAC_N_MIN
% T_PENNING


%same procesure but without for loop, not really helping...

N0=n0*ones(1,N);
[PF,eDen,rDen]=penningfraction(N0,den0);

 f=@(x)5.95*x.^5;                % This is the penning fraction distribution
 np=firstn:fix(n0/sqrt(2));      % Array of n states allowed after Penn ion
 nDen=nl*0;                      % This is the distribution of penning fraction
 ind=1:length(np); 
 nDen(ind)=f(np/n0)/sum(f(np/n0));% dividing by the sum normalizes the function
 NDEN=nDen*eDen;
 NDEN(nl==n0,:)=rDen;             % set n0 to rden
 EDEN=eDen';
 
 D_DEAC_N_MIN=NDEN(1:n_min,:);     % allow n<=n_min to decay, contributes to initial N + O from each n<=n_min on each shell     
 DEAC=sum(D_DEAC_N_MIN)';          % allow n<=n_min to decay, contributes to initial total N + O of each shell 
 NDEN(1:n_min,:)=0;   
 NL=nl*ones(1,N);               % save values for this shell in arrays
 
 
 T_PENNING=(-Ry*den0./n0^2 + Ry*rDen./n0^2 + Ry*sum(NDEN(ind,:)./nl(ind).^2)).*1./(3/2*kB.*eDen); % by energy conservation
 
 %T_penning = sum(T_PENNING.*(volume').*den0)/sum((volume').*den0); %calculate equlibrated penning temperature by weigted average
 T=5;
 
 DEAC=sum(NDEN(1:n_min,:))';       % allow n<=n_min to decay
 D_DEAC_N_MIN=NDEN(1:n_min,:);
 NDEN(1:n_min,:)=0;
 NL=nl*ones(1,N);               % save values for this shell in arrays
 DEN0=ones(ns,1)*den0;


[ni,nf,II,minn,maxn,diffsn]=buildns(nl);

index=100;

k_n_np=knnp(ni,nf,II,minn,maxn,diffsn,T);   % nl * nl matrix
kion_one=kION(nl,T); % 1*nl row
k_tbr_one=kTBR(nl(1:index),T); % 1*index row
kpd_const=kPD(nl);

k_CT=100*k_tbr_one;  % This is just an approximation of k_CT of Hydrodynamic recombination 
k_amb=1*kDR(T);        % Loss of NO+ ions to ambipolar expansion,  approximation
kpd_slow=1*kPD_slow(nl);  % predissociation rate represented by high l Rydberg only



eden=EDEN; %N*1 colum
nden=NDEN+1; 
nden=reshape(nden,N*ns,1);


D_NDEN=zeros(ns,N);
D_NDEN_OLD=zeros(ns,N);
D_EDEN=zeros(N,1);
D_DEAC=zeros(N,1);
D_EDEN_OLD=zeros(N,1);
D_DEAC_DR=zeros(N,1);
D_DEAC_PD=zeros(ns,N);
D_DEAC_N_MIN=zeros(n_min,N);


V = zeros(N,1);
D_V= zeros(N,1);
D_v= zeros(N,1)+1;

deac_dr=zeros(N,1);     %       ...           dissociative recombination 
deac_pd=zeros(N*ns,1);
deac_n_min=zeros(N*n_min,1);



if vectorize
%% start vectorizing  
    nden=reshape(nden,[ns,N]);
    deac_pd=reshape(deac_pd, [ns,N]);
    deac_n_min=reshape(deac_n_min,[n_min,N]);
    
    d_ct=zeros(numlev,N);
    d_ct(1:index,:)= k_CT.*nden(1:index,:).*eden'.^2;  % ns*N matrix with 1:index row nonzero    
    
    d_amb=k_amb*eden  % N*1 column 
    
     
    
    d_tbr=zeros(numlev,N);
    d_tbr(1:index,:)=k_tbr_one*eden'.^3;  % ns*N matrix with 1:index row nonzero
    %d_tbr=sum(d_tbr)                     % sum over all levels collapse to 1*N row
    
    d_ion=kion_one.*nden.*eden';   % ns*N matrix
    
    d_n_np=sum(k_n_np,2).*nden.*eden'; % ns*N matrix
    
    k_np_n=k_n_np';                    % ns*ns matrix
    
    k_np_n(index+1:end, :)=0;
    
    d_np_n=k_np_n*nden.*eden';      % ns*N matrix
    
    d_n_npion=sum(k_n_np(1:index,index+1:numlev),2)'*nden(1:index,:).*eden';  %  1*N  times  1*N =  1*N row
    
    D_DEAC_DR=kDR(T)*eden.^2;       % N*1 column 
    
%    D_DEAC_PD=kpd_const.*nden;      % ns*N matrix
    D_DEAC_PD=kpd_slow.*nden;      % ns*N matrix
     
    dv=D_v;                    % N*1 column  
    
 %   D_EDEN_OLD=sum(d_ion-d_tbr)'+d_n_npion';  % N*1 colum
    D_EDEN_OLD=sum(d_ion-d_tbr-d_ct)'+d_n_npion';  % N*1 colum   
  
    D_EDEN=D_EDEN_OLD-D_DEAC_DR-d_amb-eden.*dv; % N*1 column  
    
%    d_nden=d_tbr-d_ion-d_n_np+d_np_n;  % ns*N matrix
    d_nden=d_ct+d_tbr-d_ion-d_n_np+d_np_n;  % ns*N matrix
    
    D_NDEN_OLD=d_nden;   % % ns*N matrix
    
    D_DEAC=sum(d_nden(1:n_min,:))';   % N*1 column  
    
    D_DEAC_N_MIN=d_nden(1:n_min,:);   % n_min*N matrix
    
    d_nden=d_nden-D_DEAC_PD-nden.*dv'; % ns*N matrix
    
    d_nden(1:n_min,:)=0;      % ns*N matrix
    
    D_NDEN=d_nden;            % ns*N matrix
    
    D_DEAC_DR=D_DEAC_DR-deac_dr.*dv;  % N*1 column
    
    D_DEAC_PD=D_DEAC_PD-deac_pd.*dv';  % ns*N matrix
    
    D_DEAC_N_MIN=D_DEAC_N_MIN-deac_n_min.*dv'; % n_min*N matrix
   
else 
%% start looping 
 for z=1:N

        
        h=1+(z-1)*ns;
        k=z*ns;
        
        hh=1+(z-1)*n_min;
        kk=z*n_min;
        
        d_tbr=zeros(numlev,1);
        d_tbr(1:index)=k_tbr_one*eden(z)^3; % units [d_tbr] = um^-3 ns^-1
    
        d_ion=kion_one.*nden(h:k)*eden(z); % units [d_ion] = um^-3 ns^-1
        
        % rate for transfer from n to n', unit [kn_np] = um^3 ns^-1
        
        d_n_np=sum(k_n_np,2).*nden(h:k)*eden(z);            %[um^-3 ns^-1] 
        
        % rate for transfer from n' to n, [knp_n] = ns^-1
        k_np_n=zeros(numlev,numlev);
       % only to levels <= nc
         k_np_n(1:index,1:numlev)=(k_n_np(1:numlev,1:index).*(nden(h:h+numlev-1)))'; 
      
        d_np_n=sum(k_np_n,2)*eden(z);                  %[um^-3 ns^-1]
    
        % transfer from n's above ncrit(T) to eden
        k_n_npion=zeros(numlev,1);
        if index<=numlev 
            k_n_npion(1:index)=sum(k_n_np(1:index,index+1:numlev),2).*nden(h:h+index-1); % [kn_npion] = ns^-1
        end
        d_n_npion=sum(k_n_npion)*eden(z);              %[um^-3 ns^-1]
        
        %comment to deactivate DISSOCIATIVE RECOMBINATION
        D_DEAC_DR(z)=kDR(T)*eden(z)^2;
        
        %comment to deactivate PREDISSOCIATION
        D_DEAC_PD(:,z)=kpd_const.*nden(h:k); 
        
        %dv=(D_V(z)/V(z));
        dv=D_v(z);
        
        D_EDEN_OLD(z)=sum(d_ion-d_tbr)+d_n_npion; 
        D_EDEN(z)=D_EDEN_OLD(z)-D_DEAC_DR(z)-eden(z)*dv;

        d_nden=d_tbr-d_ion-d_n_np+d_np_n;
        D_NDEN_OLD(:,z)=d_nden; %array whithout radiative decay for temperature calculation
            % Implements radiative decay/PD for levels n <= n_min:
        D_DEAC(z)=sum(d_nden(1:n_min));     % aden is the number of Ry's decayed radiatively to the groundstate
        D_DEAC_N_MIN(:,z)=d_nden(1:n_min);

        d_nden=d_nden-D_DEAC_PD(:,z)-nden(h:k)*dv; %reduce by predissociation
        d_nden(1:n_min)=zeros(n_min,1); 
        
        D_NDEN(:,z)=d_nden;
        
        D_DEAC_DR(z)=D_DEAC_DR(z)-deac_dr(z)*dv;
        D_DEAC_PD(:,z)=D_DEAC_PD(:,z)-deac_pd(h:k)*dv;
        D_DEAC_N_MIN(:,z)=D_DEAC_N_MIN(:,z)-deac_n_min(hh:kk)*dv;
        
 end
%%
end



