function Y = generateY(bus,line)
% Purpose:   build admittance matrix Y from the line data
%
% Input:     bus  - bus data
%            line - line data
%
% Output:    Y    - admittance matrix
    [nb,mb]=size(bus);
    [nl,ml]=size(line);

    Y=zeros(nb,nb);         % 对导纳矩阵赋初值0
    for k=1:nl
        I=line(k,1);                       %读入线路参数
        J=line(k,2);
        Zt=line(k,3)+j*line(k,4);
        Yt=1/Zt;
        Ym=line(k,5)+j*line(k,6);
        K=line(k,7);

        if (K==0)&&(J~=0)                 % 普通线路: K=0,J~=0;
            Y(I,I)=Y(I,I)+Yt+Ym;
            Y(J,J)=Y(J,J)+Yt+Ym;
            Y(I,J)=Y(I,J)-Yt;
            Y(J,I)=Y(I,J);
        end
        if (K==0)&&(J==0)               % 对地支路: K=0,J=0,R=X=0;
            Y(I,I)=Y(I,I)+Ym;
        end
        % K<0 变压器线路: Zt和Ym为折算到K侧的值,K在i侧
        if K<0
            K=-1/K;
        end
        % K>0 变压器支路: Zt和Ym为折算到i侧的值,K在j侧
        if (K~=0)
            Y(I,I)=Y(I,I)+Yt+Ym;
            Y(J,J)=Y(J,J)+Yt/K/K;
            Y(I,J)=Y(I,J)-Yt/K;
            Y(J,I)=Y(I,J);
        end
    end
    
    
end

