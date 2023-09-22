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

load('../data/chflow/chflow_response_projection.mat')

error_interp_U = zeros(length(rs),1);   % relative error defined in terms of the squared Frobenius norm
error_interp_Psi = zeros(length(rs),1);
error_interp_rand_U = zeros(length(rs),1);
error_interp_rand_Psi = zeros(length(rs),1);

sensors_U = cell(length(rs),1);
sensors_Psi = cell(length(rs),1);
sensors_random = cell(length(rs),1);

Re_stresses_interp_U = zeros(length(rs),Ny,6);
Re_stresses_interp_Psi = zeros(length(rs),Ny,6);

tic;

for k = 1:length(rs)
    tic;
    r = rs(k);

    Ur = U(:,1:r);
    [~,sensors] = QR_sensors(Ur,r);
    not_sensed = setdiff(1:n,sensors);
    rand_sensors = not_sensed(randperm(n-r,2*r));
    sensors = [sensors rand_sensors];
    X_interp_U = Ur*(Ur(sensors,:)\X(sensors,:));
    error_interp_U(k) = norm(X-X_interp_U,'fro')^2;
    Re_stresses_interp_U(k,:,:) = chflow_Re_stresses(Fi.*X_interp_U,u_tau,Nx,Ny,Nz);
    sensors_U{k} = sensors;

    Psir = Psi(:,1:r);
    [~,sensors] = QR_sensors(Psir,r);
    not_sensed = setdiff(1:n,sensors);
    rand_sensors = not_sensed(randperm(n-r,2*r));
    sensors = [sensors rand_sensors];
    X_interp_Psi = Psir*(Psir(sensors,:)\X(sensors,:));
    error_interp_Psi(k) = norm(X-X_interp_Psi,'fro')^2;
    Re_stresses_interp_Psi(k,:,:) = chflow_Re_stresses(Fi.*X_interp_Psi,u_tau,Nx,Ny,Nz);
    sensors_Psi{k} = sensors;

    sensors = randperm(n,3*r);
    X_interp_rand_U = Ur*(Ur(sensors,:)\X(sensors,:));
    error_interp_rand_U(k) = norm(X-X_interp_rand_U,'fro')^2;
    X_interp_rand_Psi = Psir*(Psir(sensors,:)\X(sensors,:));
    error_interp_rand_Psi(k) = norm(X-X_interp_rand_Psi,'fro')^2;
    sensors_random{k} = sensors;

    toc
end
toc

error_interp_U = error_interp_U./norm2X;
error_interp_Psi = error_interp_Psi./norm2X;
error_interp_rand_U = error_interp_rand_U./norm2X;
error_interp_rand_Psi = error_interp_rand_Psi./norm2X;

save('../data/chflow/chflow_interpolatory_response_projection.mat','rs','Re_stresses_interp_U','Re_stresses_interp_Psi', ...
    'error_interp_U','error_interp_Psi','error_interp_rand_U','error_interp_rand_Psi','sensors_U','sensors_Psi','sensors_random','y','-v7.3');

%% Plot interpolatory projection error
load('../data/chflow/chflow_response_projection.mat')
load('../data/chflow/chflow_interpolatory_response_projection.mat')

c1 = [0.4940 0.1840 0.5560];
c2 = [0.3660 0.5740 0.0880];
c3 = [0.85 0.6 0.1];
lw = 1.8;
aspect =1.1;
len = 400;
f1 = figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
set(f1,'Position',[-1800 1000 1.05*len len/aspect])
plot(rs,errorU,'-','color',c2,'linewidth',lw)
hold on
plot(rs,error_interp_U,'--','color',c2,'linewidth',lw)
plot(rs,error_interp_rand_U,'-.','color',c2,'linewidth',lw)
plot(rs,errorPsi,'-','color',c1,'linewidth',lw)
plot(rs,error_interp_Psi,'--','color',c1,'linewidth',lw)
plot(rs,error_interp_rand_Psi,'-.','color',c1,'linewidth',lw)
axis([20,800,0,1.6]);
set(gca,'Fontsize',22)
pbaspect([aspect 1 1])
hold off
print(gcf,'../plots/chflow_interpolatory_projection_error','-depsc')

% f1 = figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
% set(f1,'Position',[-1800 1000 1.05*len len/aspect])
% plot(100+0*rs,errorU,'-','color',c2,'linewidth',lw)
% hold on
% plot(3*rs/n*100,error_interp_U,'--','color',c2,'linewidth',lw)
% plot(3*rs/n*100,error_interp_rand_U,'-.','color',c2,'linewidth',lw)
% plot(100+0*rs,errorPsi,'-','color',c1,'linewidth',lw)
% plot(3*rs/n*100,error_interp_Psi,'--','color',c1,'linewidth',lw)
% plot(3*rs/n*100,error_interp_rand_Psi,'-.','color',c1,'linewidth',lw)
% axis([3/n*20*100,3/n*800*100,0,1.6]);
% set(gca,'Fontsize',22)
% pbaspect([aspect 1 1])
% hold off
% print(gcf,'../plots/chflow_interpolatory_projection_error_p','-depsc')

%% Plot sensor locations

c1 = [0.4940 0.1840 0.5560];
c2 = [0.3660 0.5740 0.0880];

for k = [1,5,10,40]
    plot_chflow_sensors(x,y,z,sensors_Psi{k}(1:rs(k)),'k',c1);
    set(gcf,'Renderer','painters')
    print(gcf,['../plots/chflow_Psi_sensor_locations_r=' num2str(rs(k))],'-depsc')

    plot_chflow_sensors(x,y,z,sensors_U{k}(1:rs(k)),'k',c2);
    set(gcf,'Renderer','painters')
    print(gcf,['../plots/chflow_U_sensor_locations_r=' num2str(rs(k))],'-depsc')
end

%% Plot sensor locations in boundary layer

lc = 'k';
c1 = [0.4940 0.1840 0.5560];
c2 = [0.3660 0.5740 0.0880];
mc = [lc,'m','c'];
aspect =1.25;
len = 120;
lw = 0.75;
mt = ['o','^','s'];
ys = reshape(permute(repmat(y,1,Nx,Nz,3),[2,1,3,4]),n,1);
ds = reshape(permute(repmat([1;2;3],1,Nx,Ny,Nz),[2,3,4,1]),n,1);


k=40;
sensors = sensors_Psi{k}(1:rs(k));
length(sensors)
figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
set(gcf,'Position',[-1800 1000 1.05*len len/aspect])
[~,edges] = histcounts(log10((-abs(ys(sensors))+1)*Re_tau),12);
for j=[1,3,2]
sj = sensors(ds(sensors)==j);
histogram((-abs(ys(sj))+1)*Re_tau,10.^(edges),'Normalization','count', ...
    'FaceColor',c1,'EdgeColor',mc(j),'FaceAlpha',1,'linewidth',lw)
hold on
end
set(gca, 'xscale','log')
axis([1,Re_tau,0,250])
xticks([1,10,100])
hold off
print(gcf,'../plots/chflow_Psi_sensor_location_distribution','-depsc')

aspect =2.5;
k=40;
sensors = sensors_U{k}(1:rs(k));
figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
set(gcf,'Position',[-1800 1000 1.05*len len/aspect])
[~,edges] = histcounts(log10((-abs(ys(sensors))+1)*Re_tau),12);
for j=[1,3,2]
sj = sensors(ds(sensors)==j);
histogram((-abs(ys(sj))+1)*Re_tau,10.^(edges),'Normalization','count', ...
    'FaceColor',c2,'EdgeColor',mc(j),'FaceAlpha',1,'linewidth',lw)
hold on
end
set(gca, 'xscale','log')
axis([1,Re_tau,0,125])
xticks([1,10,100])
yticks([0,100])
hold off
print(gcf,'../plots/chflow_U_sensor_location_distribution','-depsc')

aspect =1.8;
figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
set(gcf,'Position',[-1800 1000 1.05*len len/aspect])
semilogx((y+1)*Re_tau,u_mean/u_tau,'k-','linewidth',lw)
axis([1,Re_tau,0,Inf])
xticks([1,10,100])
% yticks([0,8,16])
print(gcf,'../plots/chflow_mean_profile_wall_units','-depsc')


% for k=40%[1,5,10,40]
% sensors = sensors_Psi{k}(1:rs(k));
% figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
% set(gcf,'Position',[-1800 1000 1.05*len len/aspect])
% semilogx((y+1)*Re_tau,u_mean/u_tau,'k-','linewidth',lw)
% axis([0.5,Re_tau,0,20])
% hold on
% xticks([1,10,100])
% for j=1:3
%     sj = sensors(ds(sensors)==j);
%     [~,yj] = ismember(unique(ys(sj)),y);
%     semilogx((y(abs(yj))+1)*Re_tau,u_mean(abs(yj))/u_tau,mt(j),'MarkerSize',ms,'MarkerFaceColor',c1,'MarkerEdgeColor',mc(j))
% end
% 
% sensors = sensors_U{k}(1:rs(k));
% figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
% set(gcf,'Position',[-1800 1000 1.05*len len/aspect])
% semilogx((y+1)*Re_tau,u_mean/u_tau,'k-','linewidth',lw)
% axis([0.5,Re_tau,0,20])
% hold on
% xticks([1,10,100])
% for j=1:3
%     sj = sensors(ds(sensors)==j);
%     [~,yj] = ismember(unique(ys(sj)),y);
%     semilogx((y(abs(yj))+1)*Re_tau,u_mean(abs(yj))/u_tau,mt(j),'MarkerSize',ms,'MarkerFaceColor',c2,'MarkerEdgeColor',mc(j))
% end
% end

