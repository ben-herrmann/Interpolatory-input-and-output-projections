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

%% forcing projections
X = lsim(sys,U,t)';

r = 10;
Vr = V(:,1:r);

% orthogonal projector
U1 = Vr*(Vr'*U);

% interpolatory projector
[P,sensors] = QR_sensors(Vr,r);
U2 = Vr*(Vr(sensors,:)\U(sensors,:));

B = 0*P; s = 0.6;
for i=1:r
    xa = P(:,i)'*x;
    B(:,i) = exp(-((x-xa).^2)/(2*s^2))';
end
B = F*B;
U3 = B*((B'*B)\B')*U;

% compute responses
X1 = lsim(sys,U1,t)';
X2 = lsim(sys,U2,t)';
X3 = lsim(sys,U3,t)';

%% plot forcing and response
lc = 'w'; bc = 'k';

aspect = 1.2/1.33;
len = 450;
tf = 100;
f1 = figure(1);
set(f1,'DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex', ...
    'Position',[-1200 1000 len len/aspect])
pbaspect([aspect 1 1])

subplot(4,2,1)
s = pcolor(t,x,real(Fi*U));
s.FaceColor = 'interp';
s.EdgeColor = 'none';
colormap(redblue_b(500))
Umax = max(abs(Fi*U),[],'all');
caxis(0.3*[-Umax Umax])
axis([0,tf,-35,35])
xlabel('$t$')
ylabel('$x$')
set(gca,'layer','top', 'Color',bc, 'XColor',lc, 'YColor',lc)

subplot(4,2,3)
s = pcolor(t,x,real(Fi*U1));
s.FaceColor = 'interp';
s.EdgeColor = 'none';
% colormap(redblue(500))
caxis(0.3*[-Umax Umax])
axis([0,tf,-35,35])
xlabel('$t$')
ylabel('$x$')
set(gca,'layer','top', 'Color',bc, 'XColor',lc, 'YColor',lc)

subplot(4,2,5)
s = pcolor(t,x,real(Fi*U2));
s.FaceColor = 'interp';
s.EdgeColor = 'none';
% colormap(redblue(500))
caxis(0.3*[-Umax Umax])
axis([0,tf,-35,35])
xlabel('$t$')
ylabel('$x$')
set(gca,'layer','top', 'Color',bc, 'XColor',lc, 'YColor',lc)

subplot(4,2,7)
s = pcolor(t,x,real(Fi*U3));
s.FaceColor = 'interp';
s.EdgeColor = 'none';
% colormap(redblue(500))
caxis(0.3*[-Umax Umax])
axis([0,tf,-35,35])
xlabel('$t$')
ylabel('$x$')
set(gca,'layer','top', 'Color',bc, 'XColor',lc, 'YColor',lc)

subplot(4,2,2)
s = pcolor(t,x,real(Fi*X));
s.FaceColor = 'interp';
s.EdgeColor = 'none';
% colormap(redblue(500))
Xmax = max(abs(Fi*X),[],'all');
caxis(0.3*[-Xmax Xmax])
axis([0,tf,-35,35])
xlabel('$t$')
ylabel('$x$')
set(gca,'layer','top', 'Color',bc, 'XColor',lc, 'YColor',lc)

subplot(4,2,4)
s = pcolor(t,x,real(Fi*(X-X1)));
s.FaceColor = 'interp';
s.EdgeColor = 'none';
% colormap(redblue(500))
caxis(0.3*[-Xmax Xmax])
axis([0,tf,-35,35])
xlabel('$t$')
ylabel('$x$')
set(gca,'layer','top', 'Color',bc, 'XColor',lc, 'YColor',lc)

subplot(4,2,6)
s = pcolor(t,x,real(Fi*(X-X2)));
s.FaceColor = 'interp';
s.EdgeColor = 'none';
% colormap(redblue(500))
caxis(0.3*[-Xmax Xmax])
axis([0,tf,-35,35])
xlabel('$t$')
ylabel('$x$')
set(gca,'layer','top', 'Color',bc, 'XColor',lc, 'YColor',lc)

subplot(4,2,8)
s = pcolor(t,x,real(Fi*(X-X3)));
s.FaceColor = 'interp';
s.EdgeColor = 'none';
% colormap(redblue(500))
caxis(0.3*[-Xmax Xmax])
axis([0,tf,-35,35])
xlabel('$t$')
ylabel('$x$')
set(gca,'layer','top', 'Color',bc, 'XColor',lc, 'YColor',lc)

set(f1,'Renderer','painters','Color',bc,'InvertHardcopy','off')
% print(f1,'../plots/cgl_feedforward_b','-depsc')

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
exportgraphics(f1,'../plots/cgl_feedforward_modes.eps','BackgroundColor','none');
% print(f1,'../plots/cgl_interpolatory_input_modes','-depsc');

