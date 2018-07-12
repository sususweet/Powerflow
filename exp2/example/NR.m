clear
format long
max1=100; 
eps1=1.0e-10;
eps2=1.0e-10;
x0=[1.5  1.0]';     %初值条件
x=x0;
for i=1:max1
    A=Fd(x);
    b=F(x);
    dx=A\b;           %矩阵左除，即b除以A
    x=x-dx;
    sprintf('第%d次迭代结果：\n',i)
    [x  dx  F(x)]        %在屏幕上输出每次的x(i),dx(i),F(x(i))
    if (max(abs(dx))<eps1)&(max(abs(F(x)))<eps2)
        break
    end
end 
