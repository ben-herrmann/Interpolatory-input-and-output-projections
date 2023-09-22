clear variables; close all; clc
addpath('../src'); % Add the source files to the path

% load A matrix and integration weights
load('../data/cgl_towne_subcritical_A.mat')
n = size(A,1);
I = eye(n);
% build system
sys = ss(A,I,I,0);

%% compute obs Gramian eigenvectors
Wo = gram(sys,'o');
[V,Lambda] = eig(Wo);
lambda = diag(Lambda);
[~,p] = sort(-real(lambda));
V = V(:,p); lambda = lambda(p);

%% disturbances
dt = 0.5;
m = 10000;
t = (0:m-1)*dt;

[tt,xx] = meshgrid(t,x);
[ttc,xxc] = meshgrid(linspace(t(1),t(m),m),linspace(x(1),x(n),80));
[nc,mc] = size(ttc);
Uc = randn(nc,mc).*exp(1i*2*pi*rand(nc,mc));
U = F*interp2(ttc,xxc,Uc,tt,xx,'spline');
% for j=1:m
%     U(:,j) = U(:,j)/norm(U(:,j));
% end
U = U/sqrt(dt*norm(U,'fro')^2);

%% compute prediction error

X = lsim(sys,U,t)'; % full response
norm2X = norm(X,'fro')^2;

rs = 1:1:73;

error_Phi = zeros(length(rs),1);   % relative error defined in terms of the squared Frobenius norm
error_interp_Phi = zeros(length(rs),1);
error_interp_rand_Phi = zeros(length(rs),1);

sensors_Phi = cell(length(rs),1);
sensors_random = randperm(n,3*rs(end));

for k=1:length(rs)
    tic;

    r = rs(k);
    Vr = V(:,1:r);

    % full state sensors (orthogonal projector)
    X_Phi = lsim(sys,Vr*(Vr'*U),t)';
    error_Phi(k) = norm(X-X_Phi,'fro')^2;

    % r tailored sensors
    [~,sensors] = QR_sensors(Vr,r);
%     not_sensed = setdiff(1:n,sensors);
%     rand_sensors = not_sensed(randperm(n-r,2*r));
%     sensors = [sensors rand_sensors];
    X_interp_Phi = lsim(sys,Vr*(Vr(sensors,:)\U(sensors,:)),t)';
    error_interp_Phi(k) = norm(X-X_interp_Phi,'fro')^2;
    sensors_Phi{k} = sensors;

    % 3r random sensors
    sensors = sensors_random(1:3*r);
    X_interp_rand_Phi = lsim(sys,Vr*(Vr(sensors,:)\U(sensors,:)),t)';
    error_interp_rand_Phi(k) = norm(X-X_interp_rand_Phi,'fro')^2;

    toc
end

error_Phi = error_Phi./norm2X;
error_interp_Phi = error_interp_Phi./norm2X;
error_interp_rand_Phi = error_interp_rand_Phi./norm2X;

save('../data/chflow/cgl_interpolatory_forcing_projection.mat','rs', ...
    'error_Phi','error_interp_Phi','error_interp_rand_Phi','sensors_Phi','sensors_random','-v7.3');

%% Plot interpolatory projection error
load('../data/chflow/cgl_interpolatory_forcing_projection.mat')
c1 = [0.4940 0.1840 0.5560];
% c3 = [0.85 0.6 0.1];
lw = 1.8;
aspect =2.0;
len = 420;
f1 = figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
set(f1,'Position',[-1800 1000 1.05*len len/aspect])
semilogy(rs,error_Phi,'-','color',c1,'linewidth',lw)
hold on
semilogy(rs,error_interp_Phi,'--','color',c1,'linewidth',lw)
% plot(rs,error_interp_rand_Phi,'-.','color',c1,'linewidth',lw)
axis([0,50,0.001,5]);
yticks(10.^(-3:1))
set(gca,'Fontsize',22)
pbaspect([aspect 1 1])
hold off
set(gcf,'Renderer','painters')
print(gcf,'../plots/cgl_interpolatory_projection_error','-depsc')

%% Plot sensor locations
c1 = [0.4940 0.1840 0.5560];
c3 = [0.85 0.6 0.1];
ms = 20 - (0:50)/50*20;
lw = 0.7;
aspect =2.0;
len = 420;
f1 = figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
set(f1,'Position',[-1800 1000 1.05*len len/aspect])
scatter(rs(1),x(sensors_Phi{1}),ms(1),'o','MarkerFaceColor',c1,'MarkerEdgeColor','k','LineWidth',lw)
hold on
for k=2:50
scatter(rs(k),x(sensors_Phi{k}),ms(k),'o','MarkerFaceColor',c1,'MarkerEdgeColor','k','LineWidth',lw)
end
axis([0,50,-35,35]);
ylabel("$x$")
set(gca,'Box','on','Fontsize',22)
pbaspect([aspect 1 1])
hold off
set(gcf,'Renderer','painters')
print(gcf,'../plots/cgl_interpolatory_projection_sensors','-depsc')

%% forcing projections
X = lsim(sys,U,t)';

r = 8;
Vr = V(:,1:r);

% orthogonal projector
U1 = Vr*(Vr'*U);

% interpolatory projector
[~,sensors] = QR_sensors(Vr,r);
U2 = Vr*(Vr(sensors,:)\U(sensors,:));

% compute responses
X1 = lsim(sys,U1,t)';
X2 = lsim(sys,U2,t)';

%% plot forcing and response
aspect = 1.2;
len = 450;
tf = 100;
f1 = figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
set(f1,'Position',[-1200 1000 len len/aspect])
pbaspect([aspect 1 1])

subplot(3,2,1)
s = pcolor(t,x,real(Fi*U));
s.FaceColor = 'interp';
s.EdgeColor = 'none';
% colormap default
colormap(redblue(500))
Umax = max(abs(Fi*U),[],'all');
caxis(0.3*[-Umax Umax])
axis([0,tf,-35,35])
xlabel('$t$')
ylabel('$x$')
set(gca,'layer','top')

subplot(3,2,3)
s = pcolor(t,x,real(Fi*U1));
s.FaceColor = 'interp';
s.EdgeColor = 'none';
% colormap(redblue(500))
caxis(0.3*[-Umax Umax])
axis([0,tf,-35,35])
xlabel('$t$')
ylabel('$x$')
set(gca,'layer','top')

subplot(3,2,5)
s = pcolor(t,x,real(Fi*U2));
s.FaceColor = 'interp';
s.EdgeColor = 'none';
caxis(0.3*[-Umax Umax])
axis([0,tf,-35,35])
xlabel('$t$')
ylabel('$x$')
set(gca,'layer','top')

subplot(3,2,2)
s = pcolor(t,x,real(Fi*X));
s.FaceColor = 'interp';
s.EdgeColor = 'none';
Xmax = max(abs(Fi*X),[],'all');
caxis(0.3*[-Xmax Xmax])
axis([0,tf,-35,35])
xlabel('$t$')
ylabel('$x$')
set(gca,'layer','top')

subplot(3,2,4)
s = pcolor(t,x,real(Fi*X1));
s.FaceColor = 'interp';
s.EdgeColor = 'none';
caxis(0.3*[-Xmax Xmax])
axis([0,tf,-35,35])
xlabel('$t$')
ylabel('$x$')
set(gca,'layer','top')

subplot(3,2,6)
s = pcolor(t,x,real(Fi*X2));
s.FaceColor = 'interp';
s.EdgeColor = 'none';
caxis(0.3*[-Xmax Xmax])
axis([0,tf,-35,35])
xlabel('$t$')
ylabel('$x$')
set(gca,'layer','top')

set(f1,'Renderer','painters')
print(f1,'../plots/cgl_interpolatory_input_projection','-depsc')

%% plot modes
len = 60;
aspect = 0.3;
f1 = figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
set(f1,'Position',[-1800 1000 len len/aspect])
pbaspect([aspect 1 1])

c1 = [0.4940 0.1840 0.5560];
c2 = [0.85 0.6 0.1];
ms = 5.5;
lw = 1.8; 

% subplot(1,3,1)
plot(real(Fi*Vr(:,1)),x,'-','color',c1,'linewidth',lw)
ylim([-25 25])
hold on
plot(P'*real(Fi*Vr(:,1)),P'*x,'o','color',c2,'markerfacecolor',c2,'markersize',ms)
hold off
axis off

% subplot(1,3,2)
% plot(-real(Fi*Vr(:,2)),x,'-','color',c1,'linewidth',lw)
% ylim([-25 25])
% hold on
% plot(-P'*real(Fi*Vr(:,2)),P'*x,'o','color',c2,'markerfacecolor',c2,'markersize',ms)
% hold off
% axis off
% 
% subplot(1,3,3)
% plot(-real(Fi*Vr(:,3)),x,'-','color',c1,'linewidth',lw)
% ylim([-25 25])
% hold on
% plot(-P'*real(Fi*Vr(:,3)),P'*x,'o','color',c2,'markerfacecolor',c2,'markersize',ms)
% hold off
% axis off
set(f1,'color','none')
exportgraphics(f1,'../plots/cgl_interpolatory_input_modes.eps','BackgroundColor','none');
% print(f1,'../plots/cgl_interpolatory_input_modes','-depsc');

