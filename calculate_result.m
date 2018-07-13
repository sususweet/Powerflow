function [bus_new_result, line_result] = calculate_result (Y, bus, line)
    global U;
    global theta;
    global nodenum;
    
    nPoint = length(Y);
    G=real(Y);
    B=imag(Y);

    %����ڵ㹦�ʿ�ʼ
    for  i=1:nPoint
        bus_result(i,1) = i;
        bus_result(i,2) = U(i,1);
        bus_result(i,3) = rad2deg(theta(i,1));

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
    %�ָ��ڵ��ſ�ʼ
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
    %�ָ��ڵ��Ž���

    %������·���ʿ�ʼ
    [nl,ml]=size(line);
    for k=1:nl
        I=line(k,1);                       %������·����
        J=line(k,2);
        Zt=line(k,3)+j*line(k,4);
        Yt=1/Zt;
        Ym=line(k,5)+j*line(k,6);
        K=line(k,7);

        line_result(k,1) = line(k,1);
        line_result(k,2) = line(k,2);

        sij=0;
        sji=0;

        if (K==0)&&(J~=0)                 % ��ͨ��·: K=0,J~=0;
            sij = U(I,1)^2 * (conj(Yt) + conj(Ym)) - U(I,1)*exp(j*theta(I,1)) * conj(U(J,1)*exp(j*theta(J,1))) * conj(Yt);
            sji = U(J,1)^2 * (conj(Yt) + conj(Ym)) - U(J,1)*exp(j*theta(J,1)) * conj(U(I,1)*exp(j*theta(I,1))) * conj(Yt);
        end
        if (K==0)&&(J==0)               % �Ե�֧·: K=0,J=0,R=X=0;
            sij = U(I,1)^2 * conj(Ym);
            sji = -U(I,1)^2 * conj(Ym);
        end
        % K<0 ��ѹ����·: Zt��YmΪ���㵽K���ֵ,K��i��
        if K<0
            K=-1/K;
        end
        % K>0 ��ѹ��֧·: Zt��YmΪ���㵽i���ֵ,K��j��
        if (K~=0)
            sij = U(I,1)^2 * (conj(Yt) + conj(Ym)) - U(I,1)*exp(j*theta(I,1)) * conj(U(J,1)*exp(j*theta(J,1))) * conj(Yt) / K;
            sji = U(J,1)^2 * (conj(Yt) / K / K) - U(J,1)*exp(j*theta(J,1)) * conj(U(I,1)*exp(j*theta(I,1))) * conj(Yt) / K;
        end
        line_result(k,3) = sij;
        line_result(k,4) = sji;
        line_result(k,5) = sij + sji;
    end
    %�ָ���·��ſ�ʼ
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
    %�ָ���·��Ž���
    %������·���ʽ���

end