clear
format long
% 最大迭代次数设置
max1=100; 
% 迭代精度设置
eps1=1.0e-10;
eps2=1.0e-10;
% 初值条件
x0=[1.0  1.0 1.0];    
x=x0;

% 函数关系式设置
syms x1 x2 x3
y1=3*x1 - cos(x2*x3) - 0.5;
y2=x1^2 - 81*(x2+0.1)^2 + sin(x3) +1.06;
y3=exp(-x2*x3) +20*x3 + (10 * pi - 3) / 3; 

X = [x1;x2;x3];
Y = [y1;y2;y3];
% 求雅可比矩阵
Jacob = jacobian(Y,[x1 x2 x3]);

for i=1:max1
    b=subs(Y,symvar(Y),double(x));
    A=subs(Jacob,symvar(Jacob),double(x));
    % 求偏差量，矩阵左除，即b除以A
    dx=double(A\b);                                             
	% 求新解
    x=x'-dx;
    x=x';
    sprintf('第%d次迭代结果：\n',i)
    fx = subs(Y,symvar(Y),x);
    %在屏幕上输出每次的x(i),dx(i),F(x(i))
    double([subs(X,symvar(X),x) dx fx])
    % 如果满足精度要求则退出迭代，否则继续
    if (max(abs(dx))<eps1)&&(max(abs(fx))<eps2)
        break;
    end
end

