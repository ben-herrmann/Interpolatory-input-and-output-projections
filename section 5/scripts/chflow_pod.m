clear variables; close all; clc
addpath('../src'); % Add the source files to the path

% load data
load('../data/chflow/chflow_statistics.mat');
load_dir = '../data/chflow/chflow_dns.mat';
load(load_dir,'dt','x','y','z','n');
d = [1;2;3];
Nx = length(x); Ny = length(y); Nz = length(z);
xs = reshape(repmat(x,1,Ny,Nz,3),n,1);
ys = reshape(permute(repmat(y,1,Nx,Nz,3),[2,1,3,4]),n,1);
zs = reshape(permute(repmat(z,1,Nx,Ny,3),[2,3,1,4]),n,1);
ds = reshape(permute(repmat(d,1,Nx,Ny,Nz),[2,3,4,1]),n,1);
idx = (1:n)';
i_y = zeros(n,1);
i_z = zeros(n,1);
tic;
for i=1:n
    i_y(i) = idx(xs == xs(i) & abs(ys+ys(i))<1e-5 & zs == zs(i) & ds == ds(i));
    if zs(i) == 0
        i_z(i) = i;
    else
        i_z(i) = idx(xs == xs(i) & ys == ys(i) & abs(zs+zs(i)-(z(Nz)+z(2)))<1e-5 & ds == ds(i));
    end
end
toc

F = sqrt(reshape(repmat(0.5*cheb_weights(Ny),Nx,1,Nz,3),n,1)/(Nx*Nz));
Fi = 1./F;
fileObj = matfile(load_dir);
q_mean = [reshape(permute(repmat(u_mean,1,Nx,Nz),[2,1,3]),n/3,1); zeros(n/3,1); zeros(n/3,1)];

m = 75000;
r = 1000;
dj = m/10;
P = randn(3*m,r);
Z = zeros(n,r);
for j=1:ceil(m/dj)
    tic;
    j_start = 1+(j-1)*dj;
    j_end = min(j*dj,m);
    Xj = (F.*(fileObj.X(:,j_start:j_end)-q_mean));
    Xj_y = Xj; 
    Xj_y(n/3+1:2*n/3,:) = -Xj_y(n/3+1:2*n/3,:);
    Xj_y = Xj_y(i_y,:);
    Xj_z = Xj;
    Xj_z(2*n/3+1:n,:) = -Xj_z(2*n/3+1:n,:);
    Xj_z = Xj_z(i_z,:);
    Z = Z + Xj*P(j_start:j_end,:) + Xj_y*P(j_start+m:j_end+m,:) + Xj_z*P(j_start+2*m:j_end+2*m,:);
    disp([num2str(j) ' of' num2str(ceil(m/dj))]);
    toc
end
tic;
[Q,R] = qr(Z,0);
Y = zeros(r,3*m);
toc
for j=1:ceil(m/dj)
    tic;
    j_start = 1+(j-1)*dj;
    j_end = min(j*dj,m);
    Xj = (F.*(fileObj.X(:,j_start:j_end)-q_mean));
    Xj_y = Xj; 
    Xj_y(n/3+1:2*n/3,:) = -Xj_y(n/3+1:2*n/3,:);
    Xj_y = Xj_y(i_y,:);
    Xj_z = Xj;
    Xj_z(2*n/3+1:n,:) = -Xj_z(2*n/3+1:n,:);
    Xj_z = Xj_z(i_z,:);
    Y(:,j_start:j_end) = Q'*Xj;
    Y(:,m+j_start:m+j_end) = Q'*Xj_y;
    Y(:,2*m+j_start:2*m+j_end) = Q'*Xj_z;
    toc
end
tic;
[UY,S,V] = svd(Y,'econ');
s = diag(S);
U = Q*UY;
toc

save('../data/chflow/chflow_pod.mat','U','s','V','-v7.3');
