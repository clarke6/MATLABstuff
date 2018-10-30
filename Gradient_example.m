[x,y] = meshgrid([-10:1:10],[-10:1:10]);
z = -(x.^2 + y.^2);
figure
pcolor(x,y,z);shading('flat');colorbar
set(gca,'Fontsize',14)
title('Sclar Field Example')
xlabel('X-distance')
ylabel('Y-distance')
%note, you can also use surf if you want a 3D looking plot
figure
surf(x,y,z);
set(gca,'Fontsize',14)
title('Scalar Field Example')
xlabel('X-distance')
ylabel('Y-distance')
[px,py] = gradient(z,1,1);
figure
quiver(x,y,px,py)
set(gca,'Fontsize',14)
title('Gradient Example')
xlabel('X-distance')
ylabel('Y-distance')