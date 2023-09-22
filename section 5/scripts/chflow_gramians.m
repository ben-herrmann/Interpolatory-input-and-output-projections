clear variables; close all; clc
addpath('../src'); % Add the source files to the path

load('../data/chflow/chflow_dns.mat','x','y','z');
Lx = x(2)+x(end); Lz = z(2)+z(end);
ny = length(y);
nx = 32; nz = 16;
kx = fftshift(2*pi/Lx*[-nx/2:nx/2-1]);
kz = fftshift(2*pi/Lz*[-nz/2:nz/2-1]);
load('../data/chflow/chflow_statistics.mat');
Re_c = double(Re_c);
I = eye(2*ny-4);

lc = zeros(2*ny-4,nx,nz);
lo = zeros(2*ny-4,nx,nz);
Vc = zeros(3*ny-6,2*ny-4,nx,nz);
Vo = zeros(3*ny-6,2*ny-4,nx,nz);
Kx = zeros(2*ny-4,nx,nz);
Kz = zeros(2*ny-4,nx,nz);

for i=1:nx
    for j=1:nz
        if kx(i)~=0 || kz(j)~=0
            alpha = kx(i);
            beta = kz(j);
            k = sqrt(alpha^2+beta^2);
            [A,B,C,Q] = os_sq(u_mean,alpha,beta,Re_c);
            [F,Fi] = sqrth(Q);
            sys = ss(F*A*Fi,I,I,0);

            Wc = gram(sys,'c');
            [Vc_ij,Lc] = eig(Wc);
            lc_ij = diag(Lc);
            [~,p] = sort(-real(lc_ij));
            Vc_ij = Vc_ij(:,p);
            lc(:,i,j) = lc_ij(p);

            Wo = gram(sys,'o');
            [Vo_ij,Lo] = eig(Wo);
            lo_ij = diag(Lo);
            [~,p] = sort(-real(lo_ij));
            Vo_ij = Vo_ij(:,p);            
            lo(:,i,j) = lo_ij(p);
            for mode=1:2*ny-4
                Vc(:,mode,i,j) = C*Fi*Vc_ij(:,mode);
                Vo(:,mode,i,j) = C*Fi*Vo_ij(:,mode);
                Kx(mode,i,j) = alpha;
                Kz(mode,i,j) = beta;
            end
        end
    end
end

save('../data/chflow/chflow_gramians.mat','Re_c','x','y','z','u_mean','kx','kz', ...
    'lo','lc','Vo','Vc','Re_tau','u_tau','Re_stresses','Kx','Kz','-v7.3');

%%
w = cheb_weights(ny);
Q = diag(repmat(0.5*w(2:ny-1),1,3));
round(abs(squeeze(Vc(:,1:20,1,2))'*Q*squeeze(Vc(:,1:20,1,2))),5)
round(abs(Phi(:,1:40)'*Phi(:,1:40)),6)


%%
Lx = x(end)+x(2);
Lz = z(end)+z(2);
img = imag(exp(-1i*2*pi*(x/Lx+2*z'/Lz)));%-cos(2*pi*x/Lx).*cos(2*2*pi*z'/Lz)+sin(2*pi*x/Lx).*sin(2*2*pi*z'/Lz);
figure(1)
imagesc(real(img))
% figure(2)
% imagesc(imag(img))
img_hat = fft2(img);
figure(3)
imagesc(abs(fftshift(img_hat)))
(img_hat(2,3)/(nx*nz))
(img_hat(end,end-1)/(nx*nz))


%%

r_g = (2*ny-4)*nx*nz;
% Kx = reshape(permute(repmat(kx',1,nz,2*ny-4),[3,1,2]),r_g,1);
% Kz = reshape(permute(repmat(kz',1,nx,2*ny-4),[3,2,1]),r_g,1);
Kx = Kx(:); Kz = Kz(:);

lo_g = lo(:);
lc_g = lc(:);
Vo_g = reshape(Vo(:),3*ny-6,r_g);
Vc_g = reshape(Vc(:),3*ny-6,r_g);

[~,i_o] = sort(-real(lo_g));
[~,i_c] = sort(-real(lc_g));
lo_g = lo_g(i_o);
lc_g = lc_g(i_c);
Vo_g = Vo_g(:,i_o);
Vc_g = Vc_g(:,i_c);
kx_o = Kx(i_o); kz_o = Kz(i_o);
kx_c = Kx(i_c); kz_c = Kz(i_c);

imagesc(log10(fftshift(squeeze(real(lo(1,:,:))))))%
% ;squeeze(lc(1,nx/2:-1:2,:))])))

%%
r = 1000;
n = nx*ny*nz*3;
Psi = zeros(n,r);
Phi = zeros(n,r);
sc = lc_g(1:r);
so = lo_g(1:r);
F = sqrt(reshape(repmat(0.5*cheb_weights(ny),nx,1,nz,3),n,1)/(nx*nz));
for i=1:r
    [u,v,w] = chflow_3Dmode(Vc_g(:,i),kx,kz,kx_c(i),kz_c(i));
    Psi(:,i) = F.*[u(:);v(:);w(:)];% /(norm(F.*[u(:);v(:);w(:)]));
    [u,v,w] = chflow_3Dmode(Vo_g(:,i),kx,kz,kx_o(i),kz_o(i));
    Phi(:,i) = F.*[u(:);v(:);w(:)];% /(norm(F.*[u(:);v(:);w(:)]));
end
save('../data/chflow/chflow_forcing_response_modes.mat','Psi','sc','Phi','so','-v7.3');

%%
round(abs(Phi(:,1)'*Phi(:,12)),6)
for j=1:10
    q = Vo_g(:,j);
    [~,v,~] = chflow_3Dmode(q,kx,kz,kx_o(j),kz_o(j));
    figure()
    plot_chflow_mode(real(v),x,y,z,'f');
%     print(gcf,['../plots/chflow_phi_' num2str(j)],'-depsc')
%     q = Vc_g(:,j);
%     [u,~,~] = chflow_3Dmode(q,kx,kz,kx_c(j),kz_c(j));
%     figure()
%     plot_chflow_mode(real(u),x,y,z,'r');
%     print(gcf,['../plots/chflow_psi_' num2str(j)],'-depsc')
end
%%
j = 500;
[u,v,w] = chflow_3Dmode(Vo_g(:,j),kx,kz,kx_o(j),kz_o(j));
q = [u(:);v(:);w(:)];
wV = reshape(repmat(0.5*cheb_weights(ny),nx,1,nz,3),length(q),1)/(nx*nz);
q'*(wV.*q)