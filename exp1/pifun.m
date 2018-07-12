function y = pifun(n)
    tic;
    sum=0;
    for i = 1:n
        sum = sum + (-1)^(i-1) /(2*i-1);
    end
    y = sum * 4;
    t=toc;
    disp(['Time for pifun method:',num2str(t)])
