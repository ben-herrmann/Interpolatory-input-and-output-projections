clear variables; close all; clc
addpath('../src'); % Add the source files to the path


% load data
load('../data/chflow/chflow_statistics.mat');
load_dir = '../data/chflow/chflow_dns.mat';
load(load_dir,'dt','x','y','z','n');
Nx = length(x); Ny = length(y); Nz = length(z);
F = sqrt(reshape(repmat(0.5*cheb_weights(Ny),Nx,1,Nz,3),n,1)/(Nx*Nz));
Fi = 1./F;
fileObj = matfile(load_dir);
q_mean = [reshape(permute(repmat(u_mean,1,Nx,Nz),[2,1,3]),n/3,1); zeros(n/3,1); zeros(n/3,1)];
load('../data/chflow/chflow_pod.mat');
load('../data/chflow/chflow_forcing_response_modes.mat');
X = (F.*(fileObj.X(:,1:50:75000)-q_mean));

%%

rs = 20:20:800;

errorU = zeros(length(rs),1);   % relative error defined in terms of the squared Frobenius norm
errorPsi = zeros(length(rs),1);

Re_stresses_U = zeros(length(rs),Ny,6);
Re_stresses_Psi = zeros(length(rs),Ny,6);

tic;
X = (F.*(fileObj.X(:,1:50:75000)-q_mean));
norm2X = norm(X,'fro')^2;
Re_stresses_X = chflow_Re_stresses(Fi.*X,u_tau,Nx,Ny,Nz);

for k = 1:length(rs)
    tic;

    Ur = U(:,1:rs(k));
    XU = Ur*(Ur'*X);
    errorU(k) = norm(X-XU,'fro')^2;
    Re_stresses_U(k,:,:) = chflow_Re_stresses(Fi.*XU,u_tau,Nx,Ny,Nz);

    Psir = Psi(:,1:rs(k));
    XPsi = Psir*(Psir'*X);
    errorPsi(k) = norm(X-XPsi,'fro')^2;
    Re_stresses_Psi(k,:,:) = chflow_Re_stresses(Fi.*XPsi,u_tau,Nx,Ny,Nz);

    toc
end
toc

errorU = errorU./norm2X;
errorPsi = errorPsi./norm2X;

save('../data/chflow/chflow_response_projection.mat','rs','Re_stresses_X','Re_stresses_U','Re_stresses_Psi', ...
    'errorU','errorPsi','norm2X','y','u_tau','u_mean','Re_tau','Re_c','-v7.3');

%% Load data for plots
load('../data/chflow/chflow_response_projection.mat')


%% Plot projection error
c1 = [0.4940 0.1840 0.5560];
c2 = [0.3660 0.5740 0.0880];
lw = 1.0;
aspect =1.3;
len = 180;
f1 = figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
set(f1,'Position',[-1800 1000 1.05*len len/aspect])
plot(rs,errorU,'-','color',c2,'linewidth',lw)
hold on
plot(rs,errorPsi,'-','color',c1,'linewidth',lw)
axis([20,800,0,0.71]);
% xlabel('basis size, $r$')
% ylabel('$\mathcal{L}_2$-error')
yticks((0:0.2:0.7)')
xticks([20 (200:200:800)])
pbaspect([aspect 1 1])
hold off
print(gcf,'../plots/chflow_projection_error','-depsc')


%% Plot projected rms profiles
c1 = [0.4940 0.1840 0.5560];
c2 = [0.3660 0.5740 0.0880];
aspect =1.0;
len = 100;
lw = 0.75;
f1 = figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
set(f1,'Position',[-1800 1000 1.05*len len/aspect])
plot(y,sqrt(Re_stresses(:,1)),'k-','linewidth',lw)
hold on
plot(y,sqrt(Re_stresses(:,4)),'k-','linewidth',lw)
for k = [1,5,10,40]
    plot(y,sqrt(Re_stresses_Psi(k,:,1)),'-','color',c1,'linewidth',lw)
    plot(y,sqrt(Re_stresses_Psi(k,:,4)),'-','color',c1,'linewidth',lw)
    axis([-1.04,1.04,0,3.5])
end
hold off
print(f1,'../plots/chflow_Psi_projected_rms','-depsc')
%
f2 = figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
set(f2,'Position',[-1800 1000 1.05*len len/aspect])
semilogx((y+1)*Re_tau,sqrt(Re_stresses(:,1)),'k-','linewidth',lw)
hold on
semilogx((y+1)*Re_tau,sqrt(Re_stresses(:,4)),'k-','linewidth',lw)
for k = [1,5,10,40]
    semilogx((y+1)*Re_tau,sqrt(Re_stresses_Psi(k,:,1)),'-','color',c1,'linewidth',lw)
    semilogx((y+1)*Re_tau,sqrt(Re_stresses_Psi(k,:,4)),'-','color',c1,'linewidth',lw)
    axis([0.5,Re_tau,0,3.5])
end
xticks([1,10,100])
hold off
print(gcf,'../plots/chflow_Psi_projected_rms_wall','-depsc')
f3 = figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
set(f3,'Position',[-1800 1000 1.05*len len/aspect])
plot(y,sqrt(Re_stresses(:,1)),'k-','linewidth',lw)
hold on
plot(y,sqrt(Re_stresses(:,4)),'k-','linewidth',lw)
for k = [1,5,10,40]
    plot(y,sqrt(Re_stresses_U(k,:,1)),'color',c2,'linewidth',lw)
    plot(y,sqrt(Re_stresses_U(k,:,4)),'color',c2,'linewidth',lw)
    axis([-1.04,1.04,0,3.5])
end
hold off
print(gcf,'../plots/chflow_U_projected_rms','-depsc')
%
f4 = figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
set(f4,'Position',[-1800 1000 1.05*len len/aspect])
semilogx((y+1)*Re_tau,sqrt(Re_stresses(:,1)),'k-','linewidth',lw)
hold on
semilogx((y+1)*Re_tau,sqrt(Re_stresses(:,4)),'k-','linewidth',lw)
for k = [1,5,10,40]
    semilogx((y+1)*Re_tau,sqrt(Re_stresses_U(k,:,1)),'color',c2,'linewidth',lw)
    semilogx((y+1)*Re_tau,sqrt(Re_stresses_U(k,:,4)),'color',c2,'linewidth',lw)
    axis([0.5,Re_tau,0,3.5])
end
xticks([1,10,100])
hold off
print(gcf,'../plots/chflow_U_projected_rms_wall','-depsc')

%%
j = 55;
d = 1;
q = reshape(Fi.*U(:,j),Nx,Ny,Nz,3);
plot_chflow_mode(q(:,:,:,d),x,y,z,'f');

%% Plot projected snapshots
j = 10;
for k=[1,5,10,40]

plot_chflow_snapshot(x,y,z,Fi.*Psi(:,1:rs(k))*(Psi(:,1:rs(k))'*(X(:,j)))+q_mean,0,'k');
set(gcf,'Renderer','opengl')
print(gcf,['../plots/chflow_Psi_projected_snapshot_r=' num2str(rs(k))],'-depsc')
disp(norm(X(:,j)-Psi(:,1:rs(k))*(Psi(:,1:rs(k))'*(X(:,j))))^2/norm(X(:,j)));


plot_chflow_snapshot(x,y,z,Fi.*U(:,1:rs(k))*(U(:,1:rs(k))'*(X(:,j)))+q_mean,0,'k');
set(gcf,'Renderer','opengl')
print(gcf,['../plots/chflow_U_projected_snapshot_r=' num2str(rs(k))],'-depsc')
disp(norm(X(:,j)-U(:,1:rs(k))*(U(:,1:rs(k))'*(X(:,j))))^2/norm(X(:,j)));

end
%%
j = 10;
plot_chflow_snapshot(x,y,z,Fi.*X(:,j)+q_mean,u_mean,'w');
set(gcf,'InvertHardcopy','off');
print(gcf,'../plots/chflow_snapshot_u_mean_b','-depsc');
close all
plot_chflow_snapshot(x,y,z,Fi.*X(:,j)+q_mean,u_mean,'w');
colorbar('TickLabelInterpreter','Latex','Linewidth',1.0,'fontsize',14,'color','w')
set(gcf,'Position',[-2500 1000 400 400*0.6],'InvertHardcopy','off');
print(gcf,'../plots/chflow_snapshot_u_mean_cbar_b','-depsc');

%%
semilogx((y+1)*Re_tau,u_mean/u_tau)
axis([0.1,Re_tau,0,1/u_tau])
