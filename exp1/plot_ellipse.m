% 角度取值初始化
theta=0:0.01:2*pi;
% 计算角度对应的x, y坐标
x=5*cos(theta);
y=3*sin(theta);
% 画图，设置图形格式
plot(x,y,'-r');
title('Ellipse');
xlabel('x');
ylabel('y');
