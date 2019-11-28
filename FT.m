function [out] = FT2Dc(in)
[Nx Ny] = size(in);
h = zeros(Nx,Ny);
[m,h1]=meshgrid(1:1:Nx);
[h2,n]=meshgrid(1:1:Ny);
if Ny<=Nx
    for i=1:Nx-Ny
        h2(Ny+i,:)=h2(Ny,:);
    end
    h11=h1(:,1:Ny);
else
    h11=h1;
end

h = exp(1i.*pi.*(h11 + h2));
FT = fft2(h.*in);
out = h.*FT;

