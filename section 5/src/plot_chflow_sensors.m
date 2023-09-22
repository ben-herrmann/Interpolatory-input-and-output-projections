function fhandle = plot_chflow_sensors(x,y,z,sensors,lc,c)

% lc = 'k' or 'w'

ms = 100;
c3 = [0.85 0.6 0.1];

ymin = -1;
Nx = length(x); Ny = length(y); Nz = length(z); n = 3*Nx*Ny*Nz;

lw = 1.2;
plot3([x(1), x(end)],[z(end), z(end)],[ymin,ymin],lc,'LineWidth',lw)
axis([x(1),x(end),z(1),z(end),ymin,y(1)]);
pbaspect([x(end),z(end),y(1)-ymin])
xticks([]); yticks([]); zticks([])
set(gca,'Box','on','LineWidth',lw,'BoxStyle','full','Layer','top',...
    'XColor',lc,'YColor',lc,'ZColor',lc)

hold on
plot3([x(end), x(end)],[z(1), z(end)],[ymin,ymin],lc,'LineWidth',lw)
plot3([x(1), x(end)],[z(end), z(end)],[y(1),y(1)],lc,'LineWidth',lw)
plot3([x(end), x(end)],[z(1), z(end)],[y(1),y(1)],lc,'LineWidth',lw)
plot3([x(end), x(end)],[z(end), z(end)],[ymin,y(1)],lc,'LineWidth',lw)
plot3([x(end), x(end)],[z(1), z(1)],[ymin,y(1)],lc,'LineWidth',lw)

xs = reshape(repmat(x,1,Ny,Nz,3),n,1);
ys = reshape(permute(repmat(y,1,Nx,Nz,3),[2,1,3,4]),n,1);
zs = reshape(permute(repmat(z,1,Nx,Ny,3),[2,3,1,4]),n,1);
ds = reshape(permute(repmat([1;2;3],1,Nx,Ny,Nz),[2,3,4,1]),n,1);

mt = ['o','^','s'];
mc = [lc,'m','c'];
for j=1:3
    sj = sensors(ds(sensors)==j);
    scatter3(xs(sj),zs(sj),ys(sj),ms,mt(j),'MarkerFaceColor',c,'MarkerEdgeColor',mc(j),'LineWidth',lw)
end
hold off

if lc=='k'
    bc ='w';
else
    bc = 'k';
end
set(gcf,'Color',bc)
fhandle = gcf;
end