% This script is executed by screentshot function
% need matlab to be configured to launch the program

Vec_Unp1;
Un1(60,60)=zeros;
Un1(:)=Unp1;
clear Unp1;
Un1=transpose(Un1);
Vec_Vnp1;
Vn1(60,60)=zeros;
Vn1(:)=Vnp1;
clear Vnp1;
Vn1=transpose(Vn1);
Vec_Pnp1;
Pn1(60,60)=zeros;
Pn1(:)=Pnp1;
clear Pnp1;
Pn1=transpose(Pn1);
Vec_phinp1;
phi(60,60)=zeros;
phi(:)=phinp1;
clear phinp1;
phi=transpose(phi);
h1 = figure(1);
quiver(Un1,Vn1); 
hold on
contour(phi,[0,0],'r'); axis equal tight;
text(59,2,[' Time = ',num2str(time),' sec'],...
	'VerticalAlignment','bottom',...
	'HorizontalAlignment','right',...
	'FontSize',12,'color','red', 'BackgroundColor','white');
title('Re = 5, Bond No. = 1');
count = sprintf('%.4d',counter);
filename = strcat('Re5Bo1/plot_',count,'.eps');  
print(h1,'-depsc', filename);
close all; clear all;
