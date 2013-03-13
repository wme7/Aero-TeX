% Service functions:

%Globals1D.m-----------------------------------------------------
% Purpose: declare global variables

global N Nfp Np K
global r x  VX
global Dr LIFT
global nx Fx Fscale
global vmapM vmapP vmapB mapB Fmask
global vmapI vmapO mapI mapO
global rx J
global rk4a rk4b rk4c
global Nfaces EToE EToF
global V invV
global NODETOL
global MassMatrix

% Low storage Runge-Kutta coefficients
rk4a = [            0.0 ...
        -567301805773.0/1357537059087.0 ...
        -2404267990393.0/2016746695238.0 ...
        -3550918686646.0/2091501179385.0  ...
        -1275806237668.0/842570457699.0];
rk4b = [ 1432997174477.0/9575080441755.0 ...
         5161836677717.0/13612068292357.0 ...
         1720146321549.0/2090206949498.0  ...
         3134564353537.0/4481467310338.0  ...
         2277821191437.0/14882151754819.0];
rk4c = [             0.0  ...
         1432997174477.0/9575080441755.0 ...
         2526269341429.0/6820363962896.0 ...
         2006345519317.0/3224310063776.0 ...
         2802321613138.0/2924317926251.0];

%Startup1d.m------------------------------------------------------
function [Nv, VX, K, EToV] = MeshGen1D(xmin,xmax,K)

% function [Nv, VX, K, EToV] = MeshGen1D(xmin,xmax,K)
% Purpose  : Generate simple equidistant grid with K elements

Nv = K+1; 

% Generate node coordinates
VX = (1:Nv);
for i = 1:Nv
  VX(i) = (xmax-xmin)*(i-1)/(Nv-1) + xmin;
end

% read element to node connectivity
EToV = zeros(K, 2);
for k = 1:K
  EToV(k,1) = k; EToV(k,2) = k+1;
end
return		 

%--------------------------------------------------------------
% Purpose : Setup script, building operators, grid, metric and connectivity for 1D solver.     

% Definition of constants
Globals1D; NODETOL = 1e-10;
Np = N+1; Nfp = 1; Nfaces=2;

% Compute basic Legendre Gauss Lobatto grid
r = JacobiGL(0,0,N);

% Build reference element matrices
V  = Vandermonde1D(N, r); invV = inv(V);
Dr = Dmatrix1D(N, r, V);

% Create surface integral terms
LIFT = Lift1D();

% build coordinates of all the nodes
va = EToV (:,1)'; vb = EToV (:,2)';
x = ones(N+1,1)*VX(va) + 0.5*(r+1)*(VX(vb)-VX(va));

% calculate geometric factors
[rx,J] = GeometricFactors1D(x,Dr);

% Compute masks for edge nodes
fmask1 = find ( abs(r+1) < NODETOL)';
fmask2 = find ( abs(r-1) < NODETOL)';
Fmask  = [fmask1;fmask2]';
Fx = x(Fmask(:), :);

% Build surface normals and inverse metric at surface
[nx] = Normals1D();
Fscale = 1./(J(Fmask,:));

% Build connectivity matrix
[EToE, EToF] = Connect1D(EToV);

% Build connectivity maps
[vmapM, vmapP, vmapB, mapB] = BuildMaps1D;
%--------------------------------------------------------------
function [LIFT] = Lift1D

% function [LIFT] = Lift1D
% Purpose  : Compute surface integral term in DG formulation

Globals1D;
Emat = zeros(Np,Nfaces*Nfp);

% Define Emat
Emat(1,1) = 1.0; Emat(Np,2) = 1.0;

% inv(mass matrix)*\s_n (L_i,L_j)_{edge_n}
LIFT = V*(V'*Emat);
return
%--------------------------------------------------------------
function [rx,J] = GeometricFactors1D(x,Dr)

% function [rx,J] = GeometricFactors1D(x,Dr)
% Purpose  : Compute the metric elements for the local mappings of the 1D elements     

xr  = Dr*x; J = xr; rx = 1./J; 
return
%--------------------------------------------------------------
function [EToE, EToF] = Connect1D(EToV)

% function [EToE, EToF] = Connect1D(EToV)
% Purpose  : Build global connectivity arrays for 1D grid based on standard 
%	         EToV input array from grid generator

Nfaces = 2;
% Find number of elements and vertices
K = size(EToV,1); TotalFaces = Nfaces*K; Nv = K+1;

% List of local face to local vertex connections
vn = [1,2];

% Build global face to node sparse array
SpFToV = spalloc(TotalFaces, Nv, 2*TotalFaces);
sk = 1;
for k=1:K
  for face=1:Nfaces
     SpFToV( sk, EToV(k, vn(face))) = 1;
     sk = sk+1;
  end
end

% Build global face to global face sparse array
SpFToF = SpFToV*SpFToV' - speye(TotalFaces);

% Find complete face to face connections
[faces1, faces2] = find(SpFToF==1);

% Convert face global number to element and face numbers
element1 = floor( (faces1-1)/Nfaces )  + 1;
face1    =   mod( (faces1-1), Nfaces ) + 1;
element2 = floor( (faces2-1)/Nfaces )  + 1;
face2    =   mod( (faces2-1), Nfaces ) + 1;

% Rearrange into Nelements x Nfaces sized arrays
ind = sub2ind([K, Nfaces], element1, face1);
EToE      = (1:K)'*ones(1,Nfaces);
EToF      = ones(K,1)*(1:Nfaces);
EToE(ind) = element2; EToF(ind) = face2;
return
%--------------------------------------------------------------
function [vmapM, vmapP, vmapB, mapB] = BuildMaps1D

% function [vmapM, vmapP, vmapB, mapB] = BuildMaps1D
% Purpose: Connectivity and boundary tables for nodes given in the K # of elements,
% 	       each with N+1 degrees of freedom.

Globals1D;

% number volume nodes consecutively
nodeids = reshape(1:K*Np, Np, K);
vmapM   = zeros(Nfp, Nfaces, K); 
vmapP   = zeros(Nfp, Nfaces, K); 

for k1=1:K
  for f1=1:Nfaces
    % find index of face nodes with respect to volume node ordering
    vmapM(:,f1,k1) = nodeids(Fmask(:,f1), k1);
  end
end

for k1=1:K
  for f1=1:Nfaces
    % find neighbor
    k2 = EToE(k1,f1); f2 = EToF(k1,f1);
    
    % find volume node numbers of left and right nodes 
    vidM = vmapM(:,f1,k1); vidP = vmapM(:,f2,k2);
    
    x1  = x(vidM); x2  = x(vidP);
    
    % Compute distance matrix
    D = (x1 -x2 ).^2;
    if (D<NODETOL) vmapP(:,f1,k1) = vidP; end;
  end
end

vmapP = vmapP(:); vmapM = vmapM(:);

% Create list of boundary nodes
mapB = find(vmapP==vmapM); vmapB = vmapM(mapB);

% Create specific left (inflow) and right (outflow) maps
mapI = 1; mapO = K*Nfaces; vmapI = 1; vmapO = K*Np;
return
%-----------------------------------------------------------------
function [rhsu] = SBBGK_RHS1D(c,u,u_eq,tau)

% function [rhsrho, rhsrhou, rhsEner] = EulerRHS1D(rho, rhou ,Ener)
% Purpose  : Evaluate RHS flux in 1D Euler

Globals1D;

% compute maximum velocity for LF flux
lm = abs(c*ones(size(u)));

% Compute fluxes
uf = c.*u; 

% compute source term
us = -(u-u_eq)/tau;

% Compute jumps at internal faces
du  =zeros(Nfp*Nfaces,K);  du(:)  =  u(vmapM) -  u(vmapP); 
duf =zeros(Nfp*Nfaces,K); duf(:)  = uf(vmapM) - uf(vmapP);
LFc =zeros(Nfp*Nfaces,K); LFc(:)  = max(lm(vmapP),lm(vmapM));

% Compute fluxes at interfaces
duf(:) = nx (:).*duf (:)/2.0-LFc (:)/2.0.*du(:); 

% Boundary conditions 

% Dirichlet
% uin    = 0.0;   
% uout   = 0.0;   
% 
% % Set fluxes at inflow/outflow
% ufin =uin.^2/2; 
% lmI=lm (vmapI)/2; nxI=nx(mapI);
% duf (mapI)=nxI*(uf (vmapI)-ufin )/2.0-lmI*(u(vmapI) -uin);  
% 
% ufout=uout.^2/2; 
% lmO=lm (vmapO)/2; nxO=nx(mapO);
% duf (mapO)=nxO*(uf(vmapO) - ufout)/2.0-lmO*(u(vmapO)- uout);  

% Neumann
duf (mapI) = 0;
duf (mapO) = 0;

% compute right hand sides of the PDE's
rhsu  = -rx.*(Dr*uf)  + LIFT*(Fscale.*duf) + ((MassMatrix^-1)*us);
return

%-----------------------------------------------------------------
function ulimit = SlopeLimitN(u);

% function ulimit = SlopeLimitN(u);
% Purpose: Apply slopelimiter (Pi^N) to u assuming u an N'th order polynomial            

Globals1D;

% Compute cell averages
uh = invV*u; uh(2:Np,:)=0; uavg = V*uh; v = uavg(1,:);

% Apply slope limiter as needed.
ulimit = u; eps0=1.0e-8;

% find end values of each element
ue1 = u(1,:); ue2 = u(end,:);

% find cell averages
vk = v; vkm1 = [v(1),v(1:K-1)]; vkp1 = [v(2:K),v(K)]; 

% Apply reconstruction to find elements in need of limiting
ve1 = vk - minmod([(vk-ue1);vk-vkm1;vkp1-vk]);
ve2 = vk + minmod([(ue2-vk);vk-vkm1;vkp1-vk]);
ids = find(abs(ve1-ue1)>eps0 | abs(ve2-ue2)>eps0);

% Check to see if any elements require limiting
if(~isempty(ids))
  % create piecewise linear solution for limiting on specified elements
  uhl = invV*u(:,ids); uhl(3:Np,:)=0; ul = V*uhl;
  
  % apply slope limiter to selected elements
  ulimit(:,ids) = SlopeLimitLin(ul,x(:,ids),vkm1(ids),vk(ids),vkp1(ids));
end
return;
%----------------------------------------------------------------
function ulimit = SlopeLimitLin(ul,xl,vm1,v0,vp1);

% function ulimit = SlopeLimitLin(ul,xl,vm1,v0,vp1);
% Purpose: Apply slopelimited on linear function ul(Np,1) on x(Np,1)
%          (vm1,v0,vp1) are cell averages left, center, and right

Globals1D;

% Compute various geometric measures
ulimit = ul; h = xl(Np,:)-xl(1,:); 
x0 = ones(Np,1)*(xl(1,:) + h/2);

hN = ones(Np,1)*h;

% Limit function
ux = (2./hN).*(Dr*ul);

% ulimit = ones(Np,1)*v0+(xl-x0).*(ones(Np,1)*minmod([ux(1,:);(vp1-v0)./h;(v0-vm1)./h]));
ulimit = ones(Np,1)*v0+(xl-x0).*(ones(Np,1)*minmodB([ux(1,:);(vp1-v0)./h;(v0-vm1)./h],0.2,h));
return

%------------------------------------------------------------------
function mfunc = minmodB(v,M,h)

% function mfunc = minmodB(v,M,h)
% Purpose: Implement the TVB modified midmod function. v is a vector

mfunc = v(1,:);
ids = find(abs(mfunc) > M*h.^2);

if(size(ids,2)>0)
  mfunc(ids) = minmod(v(:,ids));
end
return