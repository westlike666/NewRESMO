



rr = 0:1:10;
thth = (0:.05:1)*2*pi/4;
[r th] = meshgrid(rr,thth)
r=linspace(0.5,1,length(th))'.*r
x = r.*cos(th);
y = r.*sin(th);
z = 1 + x.^2 - y.^2;
surf(x,y,z)