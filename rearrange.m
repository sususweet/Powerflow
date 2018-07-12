function [oribus, oriline] = rearrange(oribus, oriline)
% REARRANGE: 给节点重新编号，按照PQ、PV 和 SW 顺序进行
% 此处显示详细说明
    [nb,mb]=size(oribus);
    [nl,ml]=size(oriline);
    global nSW;
    global nPV;
    global nPQ;
    nSW = 0;                 % number of swing bus counter
    nPV = 0;                 % number of PV bus counter
    nPQ = 0;                 % number of PQ bus counter
    
    for i = 1:nb              % nb为总节点数
        type= oribus(i,6);
        if type == 3
            nSW = nSW + 1;     % increment swing bus counter
            SW(nSW,:)=oribus(i,:);
        elseif type == 2
            nPV = nPV +1;      % increment PV bus counter
            PV(nPV,:)=oribus(i,:);
        else
            nPQ = nPQ + 1;     % increment PQ bus counter
            PQ(nPQ,:)=oribus(i,:);
        end
    end;

    oribus=[PQ;PV;SW];
    newbus=[1:nb]';
    nodenum=[newbus oribus(:,1)];
    oribus(:,1)=newbus;

    for i=1:nl
        for j=1:2
            for k=1:nb
                if oriline(i,j)==nodenum(k,2)
                    oriline(i,j)=nodenum(k,1);
                    break
                end
            end
        end
    end
end

