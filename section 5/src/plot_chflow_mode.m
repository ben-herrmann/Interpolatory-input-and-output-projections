function fhandle = plot_chflow_mode(u,x,y,z,mode,lc)

lw = 1.2;
if lc=='k'
    bc ='w';
else
    bc = 'k';
end
c = 0.6;

red = [1 0 0];
green = [0 0.8 0];
yellow = [1 1 0];
blue = [0 0 0.7];
white = [0.9 0.9 0.9];
black = [0.1 0.1 0.1];

if mode=='f'
    color1 = 0.1*green+0.9*white;
    color2 = 0.6*green+0.4*black;
else
    color1 = 0.1*red+0.9*white;
    color2 = 0.6*red+0.4*black;
end

[xx,zz,yy] = meshgrid(x,z,y);
uu = permute(u,[3,1,2])/max(abs(u(:)));

plot3([x(end),x(end),x(1),x(1),x(end)],[z(1),z(end),z(end),z(1),z(1)],[y(end),y(end),y(end),y(end),y(end)],lc,'LineWidth',lw)
axis([x(1),x(end),z(1),z(end),y(end),y(1)]);
pbaspect([x(end),z(end),y(1)-y(end)])
xticks([]); yticks([]); zticks([])
set(gca,'Box','off','LineWidth',lw,'BoxStyle','full','Layer','top',...
    'XColor',lc,'YColor',lc,'ZColor',lc,'Color',bc)
hold on
plot3([x(end),x(end),x(1),x(1),x(end)],[z(1),z(end),z(end),z(1),z(1)],[y(1),y(1),y(1),y(1),y(1)],lc,'LineWidth',lw)
plot3([x(1),x(1),x(end),x(end),x(1)],[z(end),z(end),z(end),z(end),z(end)],[y(1),y(end),y(end),y(1),y(1)],lc,'LineWidth',lw)
plot3([x(1),x(1),x(end),x(end),x(1)],[z(1),z(1),z(1),z(1),z(1)],[y(1),y(end),y(end),y(1),y(1)],lc,'LineWidth',lw)

p = patch(isosurface(xx,zz,yy,uu,c));
isonormals(xx,zz,yy,uu,p)
set(p,'FaceColor',color1,'EdgeColor','none');
daspect([1 1 1])
% view([20,25]);
camlight(0,200)
lighting gouraud
% view([-70,25]);

p = patch(isosurface(xx,zz,yy,uu,-c));
isonormals(xx,zz,yy,uu,p)
set(p,'FaceColor',color2,'EdgeColor','none');
daspect([1 1 1])
lighting gouraud
% xticks([])
% yticks([])
% zticks([])
hold off

set(gcf,'Color',bc)
fhandle = gcf;
end