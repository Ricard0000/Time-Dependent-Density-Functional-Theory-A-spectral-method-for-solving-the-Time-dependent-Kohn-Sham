function [out]=trapezc3d(Temp,hx,hy,hz,Nx,Ny,Nz)

out=8*sum(sum(sum(Temp(2:Nx,2:Ny,2:Nz))));
out=out+4*sum(sum(Temp(2:Nx,:,1)));
out=out+4*sum(sum(Temp(2:Nx,:,Nz)));
out=out+4*sum(sum(Temp(1,:,2:Nz)+Temp(Nx,:,2:Nz)))+4*sum(sum(Temp(2:Nx,1,2:Nz)+Temp(2:Nx,Ny,2:Nz)));
out=out+2*sum(Temp(1,1,2:Nz)+Temp(1,Ny,2:Nz)+Temp(Nx,1,2:Nz)+Temp(Nx,Ny,2:Nz));
out=out+2*sum(Temp(1,2:Ny,1)+Temp(1,2:Ny,Nz)+Temp(Nx,2:Ny,1)+Temp(Nx,2:Ny,Nz));
out=out+2*sum(Temp(2:Nx,1,1)+Temp(2:Nx,1,Nz)+Temp(2:Nx,Ny,1)+Temp(2:Nx,Ny,Nz));
out=out+Temp(1,1,1)+Temp(1,1,Nz)+Temp(1,Ny,1)+Temp(1,Ny,Nz);
out=out+Temp(Nx,1,1)+Temp(Nx,1,Nz)+Temp(Nx,Ny,1)+Temp(Nx,Ny,Nz);
out=out*hx*hy*hz/8;
end

