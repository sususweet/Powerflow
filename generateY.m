function Y = generateY(bus,line)
% Purpose:   build admittance matrix Y from the line data
%
% Input:     bus  - bus data
%            line - line data
%
% Output:    Y    - admittance matrix
    [nb,mb]=size(bus);
    [nl,ml]=size(line);

    Y=zeros(nb,nb);         % �Ե��ɾ��󸳳�ֵ0
    for k=1:nl
        I=line(k,1);                       %������·����
        J=line(k,2);
        Zt=line(k,3)+j*line(k,4);
        Yt=1/Zt;
        Ym=line(k,5)+j*line(k,6);
        K=line(k,7);

        if (K==0)&&(J~=0)                 % ��ͨ��·: K=0,J~=0;
            Y(I,I)=Y(I,I)+Yt+Ym;
            Y(J,J)=Y(J,J)+Yt+Ym;
            Y(I,J)=Y(I,J)-Yt;
            Y(J,I)=Y(I,J);
        end
        if (K==0)&&(J==0)               % �Ե�֧·: K=0,J=0,R=X=0;
            Y(I,I)=Y(I,I)+Ym;
        end
        % K<0 ��ѹ����·: Zt��YmΪ���㵽K���ֵ,K��i��
        if K<0
            K=-1/K;
        end
        % K>0 ��ѹ��֧·: Zt��YmΪ���㵽i���ֵ,K��j��
        if (K~=0)
            Y(I,I)=Y(I,I)+Yt+Ym;
            Y(J,J)=Y(J,J)+Yt/K/K;
            Y(I,J)=Y(I,J)-Yt/K;
            Y(J,I)=Y(I,J);
        end
    end
    
    
end

