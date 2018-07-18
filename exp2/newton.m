clear
format long
% 最大迭代次数设置
max1=100; 
% 迭代精度设置
eps1=1.0e-10;
eps2=1.0e-10;
% 初值条件
x0=[1.5  1.0];     
x=x0;

% 函数关系式设置
syms x1 x2
y1=x1+2*x2-3;
y2=2*x1^2+x2^2-5;

X = [x1;x2];
Y = [y1;y2];
% 求雅可比矩阵
Jacob = jacobian(Y,[x1 x2]);

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

