clear
format long
max1=100; 
eps1=1.0e-10;
eps2=1.0e-10;
x0=[1.5  1.0]';     %��ֵ����
x=x0;
for i=1:max1
    A=Fd(x);
    b=F(x);
    dx=A\b;           %�����������b����A
    x=x-dx;
    sprintf('��%d�ε��������\n',i)
    [x  dx  F(x)]        %����Ļ�����ÿ�ε�x(i),dx(i),F(x(i))
    if (max(abs(dx))<eps1)&(max(abs(F(x)))<eps2)
        break
    end
end 
