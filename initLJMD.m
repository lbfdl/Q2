%create initial positions and velocities
%for atomic MD calculation
%  atoms in fcc structure to start
%  velocities scaled for a given reduced temperature
%input:
%  nc = number of fcc cells 
%     n = number of atoms = 4*nc^3
%     nc=2, n=32
%     nc=3, n=108
%  tin = initial temperature
%
%output:
%  n = number of atoms
%  [sx,sy,sz] = scaled coordinates
%  [vx,vy,vz] = scaled velocities
function[n,sx,sy,sz,vx,vy,vz]= initLJMD(nc,tin)
%  created scaled coordinates in an fcc lattice
ncell=4;
x=[0 .5 0 .5];
y=[0 .5 .5 0];
z=[0 0 .5 .5];
i1=0;
n = ncell*nc^3;
sx=zeros(n,1);
sy=zeros(n,1);
sz=zeros(n,1);
vx=zeros(n,1);
vy=zeros(n,1);
vz=zeros(n,1);
for k=1:nc
    for l = 1:nc
        for m = 1:nc
            for i = 1:ncell
                i1=i1+1;
                sx(i1) = (x(i) + k-1)/nc;
                sy(i1) = (y(i) + l-1)/nc;
                sz(i1) = (z(i) + m-1)/nc;
            end
        end
    end
end
%  pick velocities from Maxwell-Boltzmann distribution
%  for any temperature we want.
%  Then we will calculate the kinetic energy and thus
%  the temperature of these atoms and then we will
%  rescale the velocities to the correct temperature
k = 0;
px = 0;
py = 0;
pz = 0;
for i=1:n
    vx(i) = sqrt(-2*log(rand))*cos(2*pi*rand);
    vy(i) = sqrt(-2*log(rand))*cos(2*pi*rand);
    vz(i) = sqrt(-2*log(rand))*cos(2*pi*rand);
    px = px + vx(i);
    py = py + vy(i);
    pz = pz + vz(i);
end
% set net momentum to zero and calcuate K
px = px/n;
py = py/n;
pz = pz/n;
for i=1:n
    vx(i) = vx(i)-px;
    vy(i) = vy(i)-py;
    vz(i) = vz(i)-pz;
    k = k + vx(i)^2 + vy(i)^2 + vz(i)^2;
end
k = .5*k;
% kinetic energy of desired temperature (tin)
kin = 3*n*tin/2;
% rescale velocities
sc=sqrt(kin/k);
for i=1:n
    vx(i) = vx(i)*sc;
    vy(i) = vy(i)*sc;
    vz(i) = vz(i)*sc;
end
    


