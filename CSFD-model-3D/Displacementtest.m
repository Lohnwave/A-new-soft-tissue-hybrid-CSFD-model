clear all;close all;
load ANSYS.mat;
Fz=[0:0.025:0.6];
Dz=zeros(size(Fz));
Dz=-Gz(Fz);
figure(1)
plot(Dz,Fz,'b-^');
axis([0,9,0,0.9]);
%%
xlabel('Displacement(mm)','fontsize',16);ylabel('Foece(N)','fontsize',16);
set(gca,'fontsize',14,'fontname','Times');%设置刻度字体大小及字体
figure(2)
Dz2=-Gz(FA);
plot(DA,FA,'b-o',Dz2,FA,'r-^');
xlabel('Displacement(mm)');ylabel('Force(N)');
