%X方向分量形变曲面参数求解
%作者：罗族==================
%时间：2016-8-22=============
function Px=DzparameterX(Dz)
Px=zeros(7,1);
a=[5.12E-06;9.59E-07;2.50E-07;-3.43E-08;4.91E-10;-1.72E-12;0];
b=[-0.0608634;-0.0033395;8.48E-05;1.84E-05;-6.13E-07;5.11E-09;0];
for i=1:7
    Px(i,1)=a(i,1)+b(i,1)*Dz;
end
Px(7,1)=6.83401E-09;
end