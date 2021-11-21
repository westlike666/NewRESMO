%This part solves the rate equations and saves results in workspace
clear; clc;
tosave=true;
for density=[0.04]   % in unit um^-3, which is 10^12*cm^-3
    tic
    N=100;%number of shells
    t_max=10*1000; % time in ns 
    %steps=100;
    single=false; % single means fixed volume all the time, does not expand(dV=0 and is not coupled in to rate equatins); has to be true for single shell 
    vectorize=true % using vectorizing method is much faster than using for loop
    
    sigma_x=1*1000;
    sigma_y=0.55*1000;  %
    sigma_z=0.70*1000;  % Gaussian width convert from mm to um defualt 0.42
 

%     sigma_z=1*1000;%Gaussian width
%     sigma_y=0.55*1000;
%     sigma_x=0.70*1000;

    n=80; %PQN
    
   
    d_p=density; %peak density in um-3

    sigma_env=5;%consider number of sigma environments

    pos=linspace(0,sigma_env*sigma_z-0.5*sigma_env*sigma_z/(N-0.5),N);
    
    pos_x=linspace(0.5*sigma_env*sigma_x/(N-0.5),sigma_env*sigma_x,N)';
    pos_y=linspace(0.5*sigma_env*sigma_y/(N-0.5),sigma_env*sigma_y,N)';
    pos_z=linspace(0.5*sigma_env*sigma_z/(N-0.5),sigma_env*sigma_z,N)';
    
    
    d=arrayfun(@(z) d_p*exp(-(z^2)/(2*sigma_z^2)),pos); %Gaissian distribution of density

    foldername=['R_PQN_',num2str(n),'othertry'];
    mkdir(foldername);
    filename=[foldername,'\','l=',num2str(N),'_shell_sim_n=',num2str(n),'_d0=', num2str(d_p), '_sigma_', num2str(sigma_x/1000), 'mm_', ...
        num2str(sigma_y/1000),'mm',num2str(sigma_z/1000),'mm', '_tfinal',num2str(t_max),'ns_',num2str(sigma_env),'sigmaenv', '_fixvolume=', num2str(single)];
    
    %solve rate equations
    [time,nden,eden,deac_n_min,deac_dr,deac_pd,Te,rx,ry,rz,vx,vy,vz,vol,y0]=shell_rate_eqn_sim(d, pos_x, pos_y, pos_z, n, t_max, single, vectorize);
    
    %reduce file size by deleting shell specific information
    eden=EDEN(eden,vol);
    nden=NDEN(nden, vol);
    deac_n_min=NDEN(deac_n_min, vol);
    deac_pd=NDEN(deac_pd, vol);
    deac_dr=EDEN(deac_dr, vol);
%%  plot  
    figure()
    plot(time,eden, 'displayname', '$NO^+ + e^-$')
    hold on 
    plot(time,nden, 'displayname', '$NO^*$')
    plot(time,deac_n_min+deac_pd+deac_dr, 'displayname', '$N(^4S)+O(^3P)$')
    hold off
    legend('interpreter', 'latex')
    xlabel('ns')
    ylabel('total number of particles')
    
    figure()
    plot(time,Te)
    xlabel('ns')
    ylabel('Temperature in K')
    
    if tosave
    save(strcat([filename, '.mat']))   
    
    %clearvars y0
    %save(strcat([filename, '_small', '.mat']))
    end 
    
    toc
    
    
    
end

fprintf("Done everything\n");









