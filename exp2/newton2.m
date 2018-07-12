clear
format long
max1=100; 
eps1=1.0e-10;
eps2=1.0e-10;
x0=[1.5 6.5 -1.5];     %初值条件
x=x0;

syms x1 x2 x3
y1=x1 - 5*x2^2 + 7*x3^2 +12;
y2=3*x1*x2 + x1*x3 -11*x1;
y3=2*x2*x3 + 40*x1; 

X = [x1;x2;x3];
Y = [y1;y2;y3];
Jacob = jacobian(Y,[x1 x2 x3]);

for i=1:max1
    b=subs(Y,symvar(Y),double(x));
    A=subs(Jacob,symvar(Jacob),double(x));
    dx=double(A\b);                                             %矩阵左除，即b除以A
    x=x'-dx;
    x=x';
    sprintf('第%d次迭代结果：\n',i)
    fx = subs(Y,symvar(Y),x);
    double([subs(X,symvar(X),x) dx fx])        %在屏幕上输出每次的x(i),dx(i),F(x(i))
    if (max(abs(dx))<eps1)&&(max(abs(fx))<eps2)
        break;
    end
end

