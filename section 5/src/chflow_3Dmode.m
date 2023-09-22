function [u,v,w] = chflow_3Dmode(q,kx,kz,alpha,beta)

Nx = length(kx); Nz = length(kz);
Ny = length(q)/3+2;
q = [0;q(1:Ny-2);0;0;q(Ny-1:2*Ny-4);0;0;q(2*Ny-3:3*Ny-6);0];
Vhat = zeros(Nx,3*Ny,Nz);
% if (alpha==0 && beta==0) || alpha==kx(Nx/2+1) || beta==kz(Nz/2+1)
%     Vhat(kx==alpha,:,kz==beta) = q;
% else
%     Vhat(kx==-alpha,:,kz==-beta) = conj(q)/sqrt(2);
%     Vhat(kx==alpha,:,kz==beta) = q/sqrt(2);
% end
Vhat(kx==alpha,:,kz==beta) = q;
if (beta<=0 && alpha>0) || (beta<0 && alpha<=0) 
    V = sqrt(2)*real(permute(ifft2(permute(Nx*Nz*Vhat,[3,1,2])),[2,3,1]));
end
if (beta>=0 && alpha<0) || (beta>0 && alpha>=0)
    V = sqrt(2)*imag(permute(ifft2(permute(Nx*Nz*Vhat,[3,1,2])),[2,3,1]));
end
u = V(:,1:Ny,:);
v = V(:,Ny+1:2*Ny,:);
w = V(:,2*Ny+1:3*Ny,:);
end