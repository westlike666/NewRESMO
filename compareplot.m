
%% plot1
    close all
    figure(1)
    plot(time,eden, 'displayname', '$NO^+ + e^-$')
    hold on 
    plot(time,nden, 'displayname', '$NO^*$')
    plot(time,deac_n_min+deac_pd+deac_dr, 'displayname', '$N(^4S)+O(^3P)$')
    hold off
    
    xlabel('ns')
    ylabel('total number of particles')
    
    figure(2)
    plot(time,Te)
    xlabel('ns')
    ylabel('Temperature in K')

%% plot2
    figure(1)
    hold on
    plot(time,eden,'--','color','b', 'displayname', '$NO^+ + e^-$')
    hold on 
    plot(time,nden,'--','color', 'r', 'displayname', '$NO^*$')
    hold on 
    plot(time,deac_n_min+deac_pd+deac_dr,'--', 'color','y', 'displayname', '$N(^4S)+O(^3P)$')    
    legend('$NO^+ + e^-$','$NO^*$','$N(^4S)+O(^3P)$','interpreter', 'latex')
    
    figure(2)
    hold on
    plot(time,Te, '--','color', 'b')
    xlabel('ns')
    ylabel('Temperature in K')
    legend('expansion', 'static')
    hold off