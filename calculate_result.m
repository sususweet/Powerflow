function [bus_new_result, line_result] = calculate_result (Y, bus, line)
    global nPQ;
    global nPV;
    global nSW;
    global U;
    global theta;
    global nodenum;
    
    nPoint = length(Y);
    G=real(Y);
    B=imag(Y);

    %计算节点功率开始
    for  i=1:nPoint
        bus_result(i,1) = i;
        bus_result(i,2) = U(i,1);
        bus_result(i,3) = (180/pi) * (theta(i,1));

        sum = 0;
        for k = 1:nPoint
            sum = sum + U(i,1)*U(k,1) *(G(i,k)*cos(theta(i)-theta(k))+B(i,k)*sin(theta(i)-theta(k)));
        end
        bus_result(i,4) = sum;

        sum = 0;
        for k = 1:nPoint
            sum = sum + U(i,1)*U(k,1) *(G(i,k)*sin(theta(i)-theta(k))-B(i,k)*cos(theta(i)-theta(k)));
        end
        bus_result(i,5) = sum;
    end
    
    for  i=1:nPQ
        bus_result(i,6) = 1;
    end
    for i=nPQ + 1:nPQ + nPV
        bus_result(i,6) = 2;
    end
    for i=nPQ + nPV + 1:nPoint
        bus_result(i,6) = 3;
    end
    
    %恢复节点编号开始
    [nb,mb]=size(bus_result);
    for i=1:nb
        for k=1:nb
            if bus_result(i,1)==nodenum(k,1)
                bus_result(i,1)=nodenum(k,2);
                break
            end
        end
    end
    for i=1:nb
        for k=1:nb
            if  bus_result(k,1) == i
                bus_new_result(i, :) = bus_result(k, :);
                break;
            end
        end
    end
    %恢复节点编号结束

    %计算线路功率开始
    [nl,ml]=size(line);
    for k=1:nl
        I=line(k,1);                       %读入线路参数
        J=line(k,2);
        Zt=line(k,3)+j*line(k,4);
        Yt=1/Zt;
        Ym=line(k,5)+j*line(k,6);
        K=line(k,7);

        line_result(k,1) = line(k,1);
        line_result(k,2) = line(k,2);

        sij=0;
        sji=0;

        if (K==0)&&(J~=0)                 % 普通线路: K=0,J~=0;
            sij = U(I,1)^2 * (conj(Yt) + conj(Ym)) - U(I,1)*exp(j*theta(I,1)) * conj(U(J,1)*exp(j*theta(J,1))) * conj(Yt);
            sji = U(J,1)^2 * (conj(Yt) + conj(Ym)) - U(J,1)*exp(j*theta(J,1)) * conj(U(I,1)*exp(j*theta(I,1))) * conj(Yt);
        end
        if (K==0)&&(J==0)               % 对地支路: K=0,J=0,R=X=0;
            sij = U(I,1)^2 * conj(Ym);
            sji = -U(I,1)^2 * conj(Ym);
        end
        % K<0 变压器线路: Zt和Ym为折算到K侧的值,K在i侧
        if K<0
            K=-1/K;
        end
        % K>0 变压器支路: Zt和Ym为折算到i侧的值,K在j侧
        if (K~=0)
            sij = U(I,1)^2 * (conj(Yt) + conj(Ym)) - U(I,1)*exp(j*theta(I,1)) * conj(U(J,1)*exp(j*theta(J,1))) * conj(Yt) / K;
            sji = U(J,1)^2 * (conj(Yt) / K / K) - U(J,1)*exp(j*theta(J,1)) * conj(U(I,1)*exp(j*theta(I,1))) * conj(Yt) / K;
        end
        line_result(k,3) = sij;
        line_result(k,4) = sji;
        line_result(k,5) = sij + sji;
    end
    %恢复线路编号开始
    [nl,ml]=size(line_result);
    for i=1:nl
        for j_count=1:2
            for k=1:nb
                if line_result(i,j_count)==nodenum(k,1)
                    line_result(i,j_count)=nodenum(k,2);
                    break
                end
            end
        end
    end
    %恢复线路编号结束
    %计算线路功率结束

end