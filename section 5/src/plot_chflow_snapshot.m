function fhandle = plot_chflow_snapshot(x,y,z,q,u_mean,lc)

% lc = 'k' or 'w'

ymin = -0.9;
[xx,zz,yy] = meshgrid(x,z,y);
Nx = length(x); Ny = length(y); Nz = length(z); n = length(q);
uu = permute(reshape(q(1:n/3,:),Nx,Ny,Nz),[3,1,2]);

lw = 1.2;
c = [0.8500 0.3250 0.0980];
di_y = 5;
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

if length(u_mean)~=1
plot3(x(1)+u_mean/1.5,0*y+z(8),y,'Color',c,'LineWidth',lw)
plot3(x(1)+0*y,0*y+z(8),y,'--','Color',c,'LineWidth',lw)

[xa,za,ya] = meshgrid(x(1),z(8),y(1:di_y:end));
quiver3(xa,za,ya,reshape(u_mean(1:di_y:end),1,1,length(ya))/1.5,0*ya,0*ya,'off','Color',c,'LineWidth',lw)
end

s = slice(xx,zz,yy,uu,x(end),z(end),ymin);
set(s,'FaceColor','interp','EdgeColor','none');
caxis([0.4 0.815])
hold off

if lc=='k'
    bc ='w';
else
    bc = 'k';
end
set(gcf,'Color',bc)
fhandle = gcf;
end



% function fhandle = plot_Re_tau_180_snapshot(x,ym,z,U0m,Uts,k)
% 
% ymin = 0.1;
% [xx,zz,yy] = meshgrid(x,z,ym);
% uu = permute(squeeze(Uts(:,:,:,k)),[2,1,3]);
% 
% lw = 1.2;
% c = [0.8500 0.3250 0.0980];
% lc = 'w';
% di_y = 6;
% plot3([x(1), x(end)],[z(end), z(end)],[ymin,ymin],lc,'LineWidth',lw)
% axis([x(1),x(end),z(1),z(end),ymin,ym(end)]);
% pbaspect([x(end),z(end),ym(end)-ymin])
% xticks([]); yticks([]); zticks([])
% set(gca,'Box','on','LineWidth',lw,'BoxStyle','full','Layer','top',...
%     'XColor',lc,'YColor',lc,'ZColor',lc)
% 
% hold on
% plot3([x(end), x(end)],[z(1), z(end)],[ymin,ymin],lc,'LineWidth',lw)
% plot3([x(1), x(end)],[z(end), z(end)],[ym(end),ym(end)],lc,'LineWidth',lw)
% plot3([x(end), x(end)],[z(1), z(end)],[ym(end),ym(end)],lc,'LineWidth',lw)
% plot3([x(end), x(end)],[z(end), z(end)],[ymin,ym(end)],lc,'LineWidth',lw)
% 
% plot3(x(1)+squeeze(U0m)/1.5,0*ym+z(16),ym,'Color',c,'LineWidth',lw)
% plot3(x(1)+0*ym,0*ym+z(16),ym,'--','Color',c,'LineWidth',lw)
% 
% [xa,za,ya] = meshgrid(x(1),z(16),ym(1:di_y:end));
% quiver3(xa,za,ya,U0m(:,:,1:di_y:end)/1.5,0*ya,0*ya,'off','Color',c,'LineWidth',lw)
% 
% 
% s = slice(xx,zz,yy,uu,x(end),z(end),ymin);
% set(s,'FaceColor','interp','EdgeColor','none');
% caxis([0.5 1.1])
% hold off
% 
% if lc=='k'
%     bc ='w';
% else
%     bc = 'k';
% end
% set(gcf,'Color',bc)
% fhandle = gcf;
% end