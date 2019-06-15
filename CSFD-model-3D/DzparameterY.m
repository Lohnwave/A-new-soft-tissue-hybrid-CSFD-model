%Y方向分量形变曲面参数求解
%作者：罗族==================
%时间：2016-8-22=============
function Py=DzparameterY(Dz)
Py=zeros(7,1);
a=[7.42E-06;0;-1.12E-05;2.81E-06;-1.75E-07;5.16E-09;-4.16E-11];
b=[-5.26E-02;0;-2.96E-03;-2.49E-05;2.44E-05;-7.52E-07;6.27E-09];
for i=1:7
    Py(i,1)=a(i,1)+b(i,1)*Dz;
end
Py(2,1)=1.52004E-08;
end