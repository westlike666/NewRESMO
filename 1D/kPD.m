function kpd=kPD(n) 
    kpd=4.13e7./(2.*pi.*n.^5).*(0.3054+(3e-5).*(n.^2-16));
end

%0.3054 comes from (2*l+1).*[0:3]
%n^2-16 comes from the sum from l=4 to l=n-1

% to consider only high Rydberg (l>3), modified into
% kpd_high=4.13e7/(2*pi*n^5)*(3e-5)*(n^2-16); put in new funtion kPD_high