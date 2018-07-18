function y = pifun(n)
    % 启动计时器
    tic;
    % 初始化统计变量
    sum=0;
    % 循环求和
    for i = 1:n
        sum = sum + (-1)^(i-1) /(2*i-1);
    end
    % 得到最终结果
    y = sum * 4;
    % 输出程序用时
    t=toc;
    disp(['Time for pifun method:',num2str(t)])
end
