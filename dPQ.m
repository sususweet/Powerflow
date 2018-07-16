function [deltaP,deltaQ] =dPQ(Y,bus)
    global nPQ;
    global nPV;
    global nSW;
    global U;
    global theta;
    global P;
    global Q;
    
    nPoint = length(Y);
    G=real(Y);
    B=imag(Y);
    
    deltaP = zeros(nPoint-nSW,1);
    deltaQ = zeros(nPQ,1);
    
    %开始计算误差向量
    for row=1:nPoint-nSW
        sum = 0;
        for loop = 1:nPoint
            sum = sum + U(row)*U(loop)*(G(row,loop)*cos(theta(row)-theta(loop))+B(row,loop)*sin(theta(row)-theta(loop)));
        end
        deltaP(row,1)=P(row,1)-sum;
    end
    
    for row=1:nPQ
        sum = 0;
        for loop = 1:nPoint
            sum = sum + U(row)*U(loop)*(G(row,loop)*sin(theta(row)-theta(loop))-B(row,loop)*cos(theta(row)-theta(loop)));
        end
        deltaQ(row,1)=Q(row,1)-sum;
    end
    %误差向量计算完成
end