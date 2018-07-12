clear
format long
max1=100; 
eps1=1.0e-10;
eps2=1.0e-10;
x0=[1.0  1.0 1.0];     %��ֵ����
x=x0;

syms x1 x2 x3
y1=3*x1 - cos(x2*x3) - 0.5;
y2=x1^2 - 81*(x2+0.1)^2 + sin(x3) +1.06;
y3=exp(-x2*x3) +20*x3 + (10 * pi - 3) / 3; 

X = [x1;x2;x3];
Y = [y1;y2;y3];
Jacob = jacobian(Y,[x1 x2 x3]);

for i=1:max1
    b=subs(Y,symvar(Y),double(x));
    A=subs(Jacob,symvar(Jacob),double(x));
    dx=double(A\b);                                             %�����������b����A
    x=x'-dx;
    x=x';
    sprintf('��%d�ε��������\n',i)
    fx = subs(Y,symvar(Y),x);
    double([subs(X,symvar(X),x) dx fx])        %����Ļ�����ÿ�ε�x(i),dx(i),F(x(i))
    if (max(abs(dx))<eps1)&&(max(abs(fx))<eps2)
        break;
    end
end
