den0=[1,1,1];
n0=49;


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
 
 T_PENNING=(-Ry*den0./n0^2 + Ry*rDen./n0^2 + Ry*sum(NDEN(ind,:)./nl(ind).^2)).*1./(3/2*kB.*eDen); % by energy conservation
 
 DEAC=sum(NDEN(1:n_min,:))';       % allow n<=n_min to decay
 D_DEAC_N_MIN=NDEN(1:n_min,:);
 NDEN(1:n_min,:)=0;
 NL=nl*ones(1,N);               % save values for this shell in arrays
 DEN0=ones(ns,1)*den0;





