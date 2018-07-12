clear
format long
max1=100; 
eps1=1.0e-10;
eps2=1.0e-10;
x0=[1.5  1.0];     %初值条件
x=x0;

syms x1 x2
y1=x1+2*x2-3;
y2=2*x1^2+x2^2-5;

X = [x1;x2];
Y = [y1;y2];
Jacob = jacobian(Y,[x1 x2]);

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

