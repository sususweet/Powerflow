function [bus_new_result, line_result] = calculate_short_result (bus, Y, U, Ib, fault)
    global nodenum;
    global nSW;
    global nPV;
    global nPQ;
    
    nPoint = length(Y);
    G=real(Y);
    B=imag(Y);

    bus_result(:,1) = [1:nPoint]';
    bus_result(:,2) = U;
    bus_result(:,4) = bus(:,7);
    bus_result(:,3) = bus(:,6);
%     for  i=1:nPQ
%         bus_result(i,3) = 1;
%     end
%     for i=nPQ + 1:nPQ + nPV
%         bus_result(i,3) = 2;
%     end
%     for i=nPQ + nPV + 1:nPoint
%         bus_result(i,3) = 3;
%     end
%     
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

    %恢复线路编号开始
    [nl,ml]=size(Ib);
   
    for i=1:nl
        for j_count=1:2
            for k=1:nb
                if Ib(i,j_count)==nodenum(k,1)
                    Ib(i,j_count)=nodenum(k,2);
                    break
                end
            end
        end
    end
    %恢复线路编号结束
    
    [nl,ml]=size(Ib);
    faultI=fault(1,1);
    faultJ=fault(1,2);
    k = 1;
    for i =1:nl - 1
        if (faultI == Ib(i,1) && faultJ == Ib(i,2)) 
            continue;
        end
        if (faultI == Ib(i,2) && faultJ == Ib(i,1)) 
            continue;
        end
        line_result(k,:) = Ib(i,:);
        k = k+1;
    end
    line_result(k,:) = Ib(nl,:);
   % Ib = line_result;
    
    %line_result = Ib;
end