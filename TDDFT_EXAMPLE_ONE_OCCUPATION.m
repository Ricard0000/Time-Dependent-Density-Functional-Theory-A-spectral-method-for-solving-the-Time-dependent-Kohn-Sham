%This script will compute the time dependent density
%for an atom whose wave function is given in
%wf-st0001.dx
%The wave function at time t=0 is computed using OCTOPUS (A C/C++/Fortran program for TDDFT)
%the pseudo potential is: v0.dx
%the t=0 hartree potential is: vh.dx
%the exchange potential is: vxc.dx
%the density is: density.dx

%This code uses atomic units. In the paper we use epsilon
%to change to other units. I will leave it to the reader to use simple
%scalings to modify this code for other units.

clear

format long

tic;
Nx=61;%Number of grid points in x-direction
Ny=61;%Number of grid points in y-direction
Nz=61;%Number of grid points in z-direction

NTDFT=200;%Number of time-steps
DFTiter=0;%Set this greater than 0 for the corrector part of the wave function corrections
%(Not implemented here but you can consult me If you have questions. I might have older code)
%(that implements this, however, for this simulation,
%the Corrector part is not very significant for this small time step)

NKS=1;%Number of kohn-sham orbitals (In this example, it is one)
ax=-15;%Domain size in x (left)
bx=15;%Domain size in x (Right)
ay=-15;%Domain size in y (left)
by=15;%Domain size in y (Right)
az=-15;%Domain size in z (left)
bz=15;%Domain size in z (Right)

ftime=2;%Final time of the simulation

hx=(bx-ax)/(Nx-1);%Step sizes
hy=(by-ay)/(Ny-1);
hz=(bz-az)/(Nz-1);
T=(ftime/NTDFT);
ht=ftime/NTDFT;

time=zeros(NTDFT,1);%Setting up grids (of course you can use linspace)
x=zeros(Nx,1);%however, the grid setup is calculated almost instantly.
y=zeros(Ny,1);
z=zeros(Nz,1);
for I=1:NTDFT
    time(I,1)=ftime*(I-1)/NTDFT;
end
for I=1:Nx
    x(I)=ax+(bx-ax)/(Nx-1)*(I-1);
end
for I=1:Ny
    y(I)=ay+(by-ay)/(Ny-1)*(I-1);
end
for I=1:Nz
    z(I)=az+(bz-az)/(Nz-1)*(I-1);
end


occupations=zeros(NKS,1);%This sets the number of occupations = $f_j$ (in the paper)
occupations(1:1,1)=1;%Since there is only one, we set the first entry to 1

%This is how you load the important information from the .dx files
%That OCTOPUS outputs for the ground state density and wave functions
delimiterIn = ' ';
headerlinesIn = 7;
filename = 'v0.dx';
A = importdata(filename,delimiterIn,headerlinesIn);
v0T=A.data(:);
filename = 'vh.dx';
A = importdata(filename,delimiterIn,headerlinesIn);
vhT=A.data(:);
filename = 'vxc.dx';
A = importdata(filename,delimiterIn,headerlinesIn);
ALDAT=A.data(:);
filename = 'density.dx';
A = importdata(filename,delimiterIn,headerlinesIn);
densityT1=A.data(:);
filename = 'wf-st0001.dx';
A = importdata(filename,delimiterIn,headerlinesIn);
KST1(:,1)=A.data(:);


%I suggest you change these to GPU arrays if you have GPU capability
%GPU arrays do improve the speed of the code.
v0=zeros(Nx,Ny,Nz);
vh=zeros(Nx,Ny,Nz);
vxc=zeros(Nx,Ny,Nz);
density1=zeros(Nx,Ny,Nz);
density2=zeros(Nx,Ny,Nz);
KS1=zeros(Nx,Ny,Nz,NKS);
KS2=zeros(Nx,Ny,Nz,NKS);

%This is reshapeing data (It is computed almost instantly despite the nested loops)
L=1;
for I=1:Nx
    for J=1:Ny
        for K=1:Nz
             v0(I,J,K)=v0T(L);
             vh(I,J,K)=vhT(L);
             vxc(I,J,K)=ALDAT(L);
             density1(I,J,K)=densityT1(L);
             KS1(I,J,K,1)=KST1(L,1);
             L=L+1;
        end
    end
end

clearvars v0T vhT densityT1 densityT2 KST1 KST2 L

g1=zeros(Nx,Ny,Nz);
g2=zeros(Nx,Ny,Nz);
g3=zeros(Nx,Ny,Nz);


%(g1,g2,g3) is $r=(x,y,z)$ where $V(r)=\alpha\cdot r$. Take care to scale
%appropriately if you make this into a unit vector
for I=1:Nx
    for J=1:Ny
        for K=1:Nz
            g1(I,J,K)=ax+(bx-ax)/(Nx-1)*(I-1);
        end
    end
end
for I=1:Nx
    for J=1:Ny
        for K=1:Nz
            g2(I,J,K)=ay+(by-ay)/(Ny-1)*(J-1);
        end
    end
end
for I=1:Nx
    for J=1:Ny
        for K=1:Nz
            g3(I,J,K)=az+(bz-az)/(Nz-1)*(K-1);
        end
    end
end
vh1=vh;

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Grid for spectral method
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
Lx=bx-ax;               % --- Domain size along x
Ly=by-ay;               % --- Domain size along y
Lz=bz-az;               % --- Domain size along z
kx = (2 * pi / (Lx)) * [0 : ((Nx-1) / 2) (- (Nx-1) / 2) : (-1)]; % --- Wavenumbers along x
ky = (2 * pi / (Ly)) * [0 : ((Ny-1) / 2) (- (Ny-1) / 2) : (-1)]; % --- Wavenumbers along y
kz = (2 * pi / (Lz)) * [0 : ((Nz-1) / 2) (- (Nz-1) / 2) : (-1)]; % --- Wavenumbers along y

[Kx, Ky, Kz]  = meshgrid(kx, ky, kz); 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
density=zeros(Nx,Ny,Nz); %Density
 for I=1:NKS
     density(:,:,:)=occupations(I,1)*conj(KS1(:,:,:,I)).*KS1(:,:,:,I)+density(:,:,:);
 end
 
 KS2a=zeros(Nx,Ny,Nz,NKS); %Kohn-sham obital spin up (This is used since its only 1 wave function)
 KS2b=zeros(Nx,Ny,Nz,NKS); %Kohn-sham obital spin down (Not used)
 
 %You can consult me for the full-spin dependent code. I might have it
 %somewhere in my hard drive.
 
 
for K=1:NTDFT
    
    Z=trapezc3d(density,hx,hy,hz,Nx,Ny,Nz);
    ALDA11=1*myALDA3d(1*density1,Nx,Ny,Nz,hx,hy,hz);
    ALDA12=0.5*myALDA3d(2*density2,Nx,Ny,Nz,hx,hy,hz);
    ALDA1=ALDA11+ALDA12;
    vh1=poisol3d(density,Nx,Ny,Nz,ax,bx,ay,by,az,bz);
    pot1=ALDA11+vh1+v0;
        a(1)=cos(pi*time(K,1));
        a(2)=-pi*sin(pi*time(K,1));
        a(3)=-pi*pi*cos(pi*time(K,1))/2;
        a(4)=pi^3*sin(pi*time(K,1))/6;
        a(5)=pi^4*cos(pi*time(K,1))/24;
        a(6)=0;
        a(7)=0;
        a(8)=0;
        a(9)=0;
        alpha=1; %This is $alpha_x$ (NOT $\alpha$!!! in $V(r)=\alpha\cdot r$).
        %For this example, we polarize in the x-direction. Thus you see in the
        %code below that only the variable g1 is used.
        
        %Equations 3.13 is the next set of lines
        w=0;
        for Itt=1:9
            for Jtt=1:Itt
                w=0.5*alpha^2*(a(Jtt)*a(Itt-Jtt+1)*(2*(Itt-1)-(Jtt-1)+4)/(((Jtt-1)+1)*((Itt-1)-(Jtt-1)+1)*((Jtt-1)+2)))*ht^((Itt-1)+3)/((Itt-1)+3)+w;
            end
        end
        T=0;
        for Itt=1:5
            T=alpha*(a(Itt)*ht^((Itt-1)+1)/((Itt-1)+1))+T;
        end
        u=0;
        for Itt=1:5
            u=alpha*(a(Itt)*ht^((Itt-1)+2))/((Itt)*(Itt+1))+u;
        end
%    figure(111)
%    mesh(y,x,real(KS1(:,:,(Nz+1)/2,1)));drawnow;
%    This loop computes the Kohn-Sham orbitals Equation 3.16
%    You might want to try parfor here. I dont quite remember how
%    much parfar speeds up the calculation.
    for I=1:NKS
        KS2a(:,:,:,I)=KS1(:,:,:,I).*exp(-1i*ht/2*(pot1));
        u_hat=fftn(KS2a(:,:,:,I));
        u_hat=u_hat.*exp(-1i*ht/2*((Kx).^2+(Ky).^2+(Kz).^2)).*exp(-1i*Ky*u);
        KS2a(:,:,:,I)=ifftn(u_hat);
        normalizationf=exp(1i*w);
        KS2a(:,:,:,I)=normalizationf*exp(-1i*T*(g1+u)).*KS2a(:,:,:,I);
        KS2a(:,:,:,I)=KS2a(:,:,:,I).*exp(-1i*ht/2*(pot1));
    end
    density12=zeros(Nx,Ny,Nz);
    density22=zeros(Nx,Ny,Nz);
    for I=1:NKS
        density12(:,:,:)=occupations(I,1)*conj(KS2a(:,:,:,I)).*KS2a(:,:,:,I)+density12(:,:,:);
    end
    
    KS1=KS2a;
%    KS2=KS2b;
    density2=density12+density22;
    density=density2;       
    
    %Calculation of polarizations
    polarizationx(K,1)=trapezc3d(density.*g1,hx,hy,hz,Nx,Ny,Nz);
    polarizationy(K,1)=trapezc3d(density.*g2,hx,hy,hz,Nx,Ny,Nz);
    polarizationz(K,1)=trapezc3d(density.*g3,hx,hy,hz,Nx,Ny,Nz);
    
    K %This tells you when the code is going to finish
    %When K reaches NTDFT, you are done with the simulation.
end

density3da=zeros(Nx,Ny,Nz);
density3db=zeros(Nx,Ny,Nz);
for I=1:NKS
    density3da(:,:,:)=occupations(I,1)*conj(KS1(:,:,:,I)).*KS1(:,:,:,I)+density3da(:,:,:);
end
%for I=1:NKS
%    density3db(:,:,:)=occupations(I,1)*conj(KS2(:,:,:,I)).*KS2(:,:,:,I)+density3db(:,:,:);
%end

density3d=density3da;%+density3db;

t=toc;
t=floor(t);
minutes=floor(t/60);
seconds=t-minutes*60;
disp('computation time for the Gaussian beam method is')
minutes
disp('minutes and')
seconds
disp('seconds')




for I=1:Nx
    x(I)=ax+(bx-ax)/(Nx-1)*(I-1);
end

for I=1:Ny
    y(I)=ay+(by-ay)/(Ny-1)*(I-1);
end

for I=1:Nz
    z(I)=az+(bz-az)/(Nz-1)*(I-1);
end



figure(1)
surf(y,x,real(density3d(:,:,(Nz+1)/2)))
colormap(jet)
%colormap(hot)
title('z=0 cross section of $\rho(r,t_f)$ using FGA$','Interpreter','Latex','FontSize',18)
figure(2)
pcolor(y,x,density3d(:,:,(Nz+1)/2))
colormap(jet)
%colormap(hot)
title('z=0 cross section of $\rho(r,t_f)$ using FGA','Interpreter','Latex','FontSize',18)
shading interp
h=colorbar;
xlabel('$x$','Interpreter','Latex','FontSize',18)
set(gca,'FontSize',16)
ylabel('$y$','Interpreter','Latex','FontSize',18)
set(gca,'FontSize',16)

figure(3)
scatter(time,polarizationx,'r','x')
title('dipole moment using FGA','Interpreter','Latex','FontSize',18)

