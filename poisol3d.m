%3d poisol

function [out_small]=poisol3d(In,Nx,Ny,Nz,ax,bx,ay,by,az,bz,g1,g2,g3)

%Nx=32;              % --- Number of Fourier harmonics along x (should be a multiple of 2)  
%Ny=32;              % --- Number of Fourier harmonics along y (should be a multiple of 2)  
%Nz=32;
%ax=-13.6;
%bx=13.6;
%ay=-13.6;%
%by=13.6;
%az=-13.6;
%bz=13.6;
%p1=boole3d(In.*g1,ax,bx,ay,by,az,bz,Nx,Ny,Nz);
%p2=boole3d(In.*g2,ax,bx,ay,by,az,bz,Nx,Ny,Nz);
%p3=boole3d(In.*g3,ax,bx,ay,by,az,bz,Nx,Ny,Nz);

hx=(bx-ax)/(Nx-1);
hy=(by-ay)/(Ny-1);
hz=(bz-az)/(Nz-1);


%Incr=ceil(Nx/2);
Incr=0;
NNx=Nx+Incr+Incr;
NNy=Ny+Incr+Incr;
NNz=Nz+Incr+Incr;

Inc_out=zeros(NNx,NNy,NNz);
Inc_out(Incr+1:Incr+Nx,Incr+1:Ny+Incr,1+Incr:Incr+Nz)=In(:,:,:);

aax=ax-Incr*hx;
bbx=bx+Incr*hx;
aay=ay-Incr*hy;
bby=by+Incr*hy;
aaz=az-Incr*hz;
bbz=bz+Incr*hz;


Lx=bbx-aax;               % --- Domain size along x
Ly=bby-aay;               % --- Domain size along y
Lz=bbz-aaz;

% --- Wavenumbers
%kx = (2 * pi / Lx) * [0 : (Nx / 2 - 1) (- Nx / 2) : (-1)]; % --- Wavenumbers along x
%ky = (2 * pi / Ly) * [0 : (Ny / 2 - 1) (- Ny / 2) : (-1)]; % --- Wavenumbers along y
%kz = (2 * pi / Lz) * [0 : (Nz / 2 - 1) (- Nz / 2) : (-1)]; % --- Wavenumbers along y

kx = (2 * pi / Lx) * [0 : ((NNx-1) / 2) (- (NNx-1) / 2) : (-1)]; % --- Wavenumbers along x
ky = (2 * pi / Ly) * [0 : ((NNy-1) / 2) (- (NNy-1) / 2) : (-1)]; % --- Wavenumbers along y
kz = (2 * pi / Lz) * [0 : ((NNz-1) / 2) (- (NNz-1) / 2) : (-1)]; % --- Wavenumbers along y

%[Kx, Ky, Kz]  = meshgrid(kx, ky, kz); 
[Kx, Ky, Kz]  = meshgrid(ky, kx, kz); 
fHat=fftn(Inc_out);


% --- Denominator of the unknown spectrum
den             = -(Kx.^2 + Ky.^2 + Kz.^2);

den(1,1,1)      = 1; % --- Avoid division by zero at wavenumber (0, 0)

% --- Unknown determination
uHat            = ifftn(fHat./ den);

% uHat(1, 1)      = 0;            % --- Force the unknown spectrum at (0, 0) to be zero
out               = real(uHat);
out               = out-out(1,1,1);   % --- Force arbitrary constant to be zero by forcing u(1, 1) = 0
%out               = out-out((NNx-1)/2,(NNy-1)/2,(NNz-1)/2);
%Integrand=zeros(NNx,NNy,NNz);
%for I=1:NNx
%    for J=1:NNy
%        for K=1:NNz
%            x=(bbx-aax)/(NNx-1)*(I-1);
%            y=(bby-aay)/(NNy-1)*(J-1);
%            z=(bbz-aaz)/(NNz-1)*(K-1);
%            Integrand(I,J,K)=Inc_out(I,J,K)/sqrt(x^2+y^2+z^2);
%        end
%    end
%end

%Integrand=zeros(Nx,Ny,Nz);
%for I=1:Nx
%    for J=1:Ny
%        for K=1:Nz
%            x=ax+(bx-ax)/(Nx-1)*(I-1);
%            y=bx+(by-ay)/(Ny-1)*(J-1);
%            z=bz+(bz-az)/(Nz-1)*(K-1);
%            Integrand(I,J,K)=In(I,J,K)/sqrt((x)^2+(y)^2+(z)^2);
%        end
%    end
%end
%Integrand((Nx+1)/2,(Ny+1)/2,(Nz+1)/2)=In((Nx+1)/2,(Ny+1)/2,(Nz+1)/2)/(10^-3);
%Integrand(1,1,1)=In(1,1,1)/(10^-4);
%L1=trapezc3d(Integrand,hx,hy,hz,Nx,Ny,Nz)

%Multipole expansion of (1,1,1)
%Z=round(trapezc3d(In,hx,hy,hz,Nx,Ny,Nz));
%r=sqrt(aax^2+aay^2+aaz^2);
%L1=Z/r+(p1*aax+p2*aay+p3*aaz)/r^3;


out=-out*4*pi;%+L1;
out_small(:,:,:)=out(Incr+1:Incr+Nx,Incr+1:Ny+Incr,1+Incr:Incr+Nz);
out_small=out_small;
end

