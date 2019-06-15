%Z��������α�����������
%���ߣ�����==================
%ʱ�䣺2016-8-22=============
function Pz=DxparameterZ(Dx)
Pz=zeros(5,1);
a=[-8.08E-06;9.97E-07;-4.15E-08;5.72E-10;0];
b=[0.1875753;0.0090317;-0.001053;1.46E-05;0];
for i=1:5
    Pz(i,1)=a(i,1)+b(i,1)*Dx;
end
Pz(5,1)=-1.06E-09;
end