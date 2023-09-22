clear variables; close all; clc
addpath('../src'); % Add the source files to the path

% load A matrix and integration weights
load('../data/cgl_towne_subcritical_A.mat')
n = size(A,1);
I = eye(n);
% build system
sys = ss(A,I,I,0);

%% compute obs Gramian eigenvectors and balanced modes
Wo = gram(sys,'o');
[V,Lambda] = eig(Wo);
lambda = diag(Lambda);
[~,p] = sort(-real(lambda));
V = V(:,p);
Wc = gram(sys,'c');
[Vc,Lambda] = eig(Wc);
lambda = diag(Lambda);
[~,p] = sort(-real(lambda));
Vc = Vc(:,p);
[~,~,Ti,T] = balreal(sys);
Vb = Ti';
Vcb = T;
%%
j = 1;
figure(1)
plot(x,abs(Fi*V(:,j)),'-')
hold on
plot(x,abs(Fi*Vc(:,j)),'-')
plot(x,abs(Fi*Vb(:,j))/norm(Vb(:,j)),'--')
plot(x,abs(Fi*Vcb(:,j))/norm(Vcb(:,j)),'--')
xlim([-25 25])
hold off
% figure(2)
% plot(x,abs(Fi*V(:,j)).*abs(Fi*Vc(:,j)),'-')
% hold on
% plot(x,(abs(Fi*Vb(:,j))/norm(Vb(:,j))).*(abs(Fi*Vcb(:,j))/norm(Vcb(:,j))),'--')
% xlim([-25 25])
% hold off
%% disturbances and POD
dt = 0.5;
m = 10000;
t = (0:m-1)*dt;

[tt,xx] = meshgrid(t,x);
[ttc,xxc] = meshgrid(linspace(t(1),t(m),m),linspace(x(1),x(n),80));
[nc,mc] = size(ttc);
Uc = randn(nc,mc).*exp(1i*2*pi*rand(nc,mc));
U = F*interp2(ttc,xxc,Uc,tt,xx,'spline');
U = U/sqrt(dt*norm(U,'fro')^2);

[V2,~,~] = svd(U,'econ');

%% compute prediction error

X = lsim(sys,U,t)'; % full response
norm2X = norm(X,'fro')^2;

rs = 1:1:120;

error_Phi = zeros(length(rs),1);   % relative error defined in terms of the squared Frobenius norm
error_U = zeros(length(rs),1);
error_Phi_bal = zeros(length(rs),1);

for k=1:length(rs)
    tic;

    r = rs(k);
    Vr = V(:,1:r);
    V2r = V2(:,1:r);
    Vbr = Vb(:,1:r);

    % observability Gramian eigenmodes
    X_Phi = lsim(sys,Vr*(Vr'*U),t)';
    error_Phi(k) = norm(X-X_Phi,'fro')^2;

    % forcing POD
    X_U = lsim(sys,V2r*(V2r'*U),t)';
    error_U(k) = norm(X-X_U,'fro')^2;

    % adjoint balancing modes
    X_Phi_bal = lsim(sys,(Vbr/(Vbr'*Vbr))*(Vbr'*U),t)';
    error_Phi_bal(k) = norm(X-X_Phi_bal,'fro')^2;

    toc
end

error_Phi = error_Phi./norm2X;
error_U = error_U./norm2X;
error_Phi_bal = error_Phi_bal./norm2X;

save('../data/chflow/cgl_forcing_projection.mat','rs', ...
    'error_Phi','error_U','error_Phi_bal','-v7.3');

%% Plot forcing projection error
c1 = [0.4940 0.1840 0.5560];
c2 = [0.3660 0.5740 0.0880];
% c3 = [0.85 0.6 0.1];
lw = 1.8;
aspect =1.7;
len = 400;
f1 = figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
set(f1,'Position',[-1800 1000 1.05*len len/aspect])
semilogy(rs,error_Phi,'-','color',c1,'linewidth',lw)
hold on
semilogy(rs,error_U,'-','color',c2,'linewidth',lw)
plot(rs,error_Phi_bal,'--','color',c1,'linewidth',lw)
axis([0,60,0.001,5]);
yticks(10.^(-3:1))
set(gca,'Fontsize',22)
pbaspect([aspect 1 1])
hold off
set(gcf,'Renderer','painters')
% print(gcf,'../plots/cgl_forcing_projection_error','-depsc')

%% compute responses
X = lsim(sys,U,t)';

p1 = 3;
U1 = V(:,1:p1)*V(:,1:p1)'*U;
p2 = 20;
U2 = V2(:,1:p2)*V2(:,1:p2)'*U;
X1 = lsim(sys,U1,t)';
X2 = lsim(sys,U2,t)';
Ub = Vcb(:,1:p1)*Vb(:,1:p1)'*U;
% Ub = (Vb(:,1:p1)/(Vb(:,1:p1)'*Vb(:,1:p1)))*Vb(:,1:p1)'*U;
Xb = lsim(sys,Ub,t)';

%% plot forcing and response
aspect = 1.2;
len = 450;
tf = 200;
f1 = figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
set(f1,'Position',[-1200 1000 len len/aspect])
pbaspect([aspect 1 1])

subplot(3,2,1)
s = pcolor(t,x,real(Fi*U));
s.FaceColor = 'interp';
s.EdgeColor = 'none';
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
colormap(redblue(500))
caxis(0.3*[-Umax Umax])
axis([0,tf,-35,35])
xlabel('$t$')
ylabel('$x$')
set(gca,'layer','top')

subplot(3,2,5)
s = pcolor(t,x,real(Fi*U2));
s.FaceColor = 'interp';
s.EdgeColor = 'none';
colormap(redblue(500))
Umax = max(abs(Fi*U),[],'all');
caxis(0.3*[-Umax Umax])
axis([0,tf,-35,35])
xlabel('$t$')
ylabel('$x$')
set(gca,'layer','top')

subplot(3,2,2)
s = pcolor(t,x,real(Fi*X));
s.FaceColor = 'interp';
s.EdgeColor = 'none';
colormap(redblue(500))
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
colormap(redblue(500))
caxis(0.3*[-Xmax Xmax])
axis([0,tf,-35,35])
xlabel('$t$')
ylabel('$x$')
set(gca,'layer','top')

subplot(3,2,6)
s = pcolor(t,x,real(Fi*X2));
s.FaceColor = 'interp';
s.EdgeColor = 'none';
colormap(redblue(500))
caxis(0.3*[-Xmax Xmax])
axis([0,tf,-35,35])
xlabel('$t$')
ylabel('$x$')
set(gca,'layer','top')

set(f1,'Renderer','painters')
print(f1,'../plots/cgl_input_projection','-depsc')

%% plot modes
len = 180;
aspect = 0.65;
f1 = figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
set(f1,'Position',[-1800 1000 len len/aspect])
pbaspect([aspect 1 1])

c1 = [0.4940 0.1840 0.5560];
c2 = [0.3660 0.5740 0.0880];
lw = 1.0; 

subplot(2,3,1)
plot(real(Fi*V(:,1)),x,'-','color',c1,'linewidth',lw)
ylim([-35 35])
axis off

subplot(2,3,2)
plot(-real(Fi*V(:,2)),x,'-','color',c1,'linewidth',lw)
ylim([-35 35])
axis off

subplot(2,3,3)
plot(-real(Fi*V(:,3)),x,'-','color',c1,'linewidth',lw)
ylim([-35 35])
axis off

subplot(2,3,4)
plot(real(Fi*V2(:,1)),x,'-','color',c2,'linewidth',lw)
ylim([-35 35])
axis off

subplot(2,3,5)
plot(real(Fi*V2(:,2)),x,'-','color',c2,'linewidth',lw)
ylim([-35 35])
axis off

subplot(2,3,6)
plot(real(Fi*V2(:,20)),x,'-','color',c2,'linewidth',lw)
ylim([-35 35])
axis off

print(f1,'../plots/cgl_input_modes','-depsc');