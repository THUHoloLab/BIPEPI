function [p] = Propagator1(M, N, lmada, pitch, z)

p = zeros(M,N);
k = 2*pi/lmada;
lx = N*pitch;
ly = M*pitch;
[x,y] = meshgrid(linspace(-lx/2,lx/2-pitch,N),linspace(-ly/2,ly/2-pitch,M));
[fx,fy] = meshgrid(linspace(-1/2/pitch,1/2/pitch-1/lx,N),linspace(-1/2/pitch,1/2/pitch-1/ly,M));
z1 = -z;
p = exp(1i*k*z*sqrt(1-(lmada*fx).^2-(lmada*fy).^2));
