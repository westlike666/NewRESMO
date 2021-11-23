function kpd_slow=kPD_slow(n) %calculating the slow predissociating rate for only summing high l  

    kpd_slow=4.13e7./(2*pi.*n.^5).*(3e-5).*(n.^2-16); %assumming l is mixed to high l, l>3
    
    kpd_slow(1:4)=kPD(n(1:4)); % reset n=1,2,3,4 to original fast kpd. Otherwise will yield negative value
    
    %kpd_slow=0*n; This is atother approximation by setting kpd to 0.
end


%n^2-16 comes from the sum from l=4 to l=n-1

% to consider only high Rydberg (l>3).
