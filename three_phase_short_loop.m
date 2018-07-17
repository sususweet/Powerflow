function Y = three_phase_short_loop(dfile)
    % 输入初始数据：节点编号、电网参数及故障参数；

    lfile =length(dfile);
    % strip off .m
    eval(dfile(1:lfile-2));
    %
    %     [dfile,pathname]=uigetfile('*.m','Select Data File');
    %     if pathname == 0
    %         error(' you must select a valid data file')
    %     else
    %         lfile =length(dfile);
    %         % strip off .m
    %         eval(dfile(1:lfile-2));
    %     end

    global nodenum;
    global theta;
    global P;
    global Q;
    tic;

    bus_ori = bus;
    %节点重新编号开始
    [bus, line] = rearrange(bus, line);
    %节点重新编号结束
    [nb,mb]=size(bus);
    [nl,ml]=size(line);

    filename = ['./result/output_', dfile, '.dat'];
    myf=fopen(filename,'w','n','UTF-8');

    filename_graph = ['./result/output_graph_', dfile, '.dat'];
    myf_graph=fopen(filename_graph,'w','n','UTF-8');

    % 形成节点导纳矩阵
    Y= generateY(bus,line);
    max_index_bus = 0;
    max_short_current = 0;

    fprintf(myf, 'Short circuit current of each node: \n');
    fprintf(myf, 'Node I    Short circuit current effective value \n');
    fprintf(myf_graph, 'Short circuit current of each node: \n');
    fprintf(myf_graph, 'Node I    Short circuit current effective value \n');
    
    for index_bus = 1:nb
        fault = [bus_ori(index_bus,1)  0  0];

        % 对节点导纳矩阵进行修正，形成修正导纳矩阵
        % 发电机节点修正开始
        fixY = Y;
        for k=1:nb
            busType=bus(k,7);
            if (busType==1)
                fixY(bus(k,1),bus(k,1)) =  fixY(bus(k,1),bus(k,1)) + 1/j/bus(k,8);
            end
        end
        % 发电机节点修正完成

        % 故障点修正开始
        [nf,mf]= size(fault);
        for k=1:nf
            nodeI=fault(k,1);
            nodeJ=fault(k,2);
            dis=fault(k,3);

            % 节点三相短路故障处理开始
            if (nodeI==0)
                fixY(nodeJ,nodeJ) = 999999+j*999999;
                continue;
            end
            if (nodeJ==0)
                fixY(nodeI,nodeI) = 999999+j*999999;
                continue;
            end
            if (dis==0)&&(nodeI*nodeJ~=0)
                fixY(nodeI,nodeI) = 999999+j*999999;
                continue;
            end
            if (dis==1)&&(nodeI*nodeJ~=0)
                fixY(nodeJ,nodeJ) = 999999+j*999999;
                continue;
            end
            % 节点三相短路故障处理结束

            % 线路三相短路故障处理开始
            if (dis~=1)&&(dis~=0)&&(nodeI*nodeJ~=0)
                fixY(nodeI,nodeI) = fixY(nodeI,nodeI) + fixY(nodeI,nodeJ) - fixY(nodeI,nodeJ)/dis;
                fixY(nodeJ,nodeJ) = fixY(nodeJ,nodeJ) + fixY(nodeI,nodeJ) - fixY(nodeI,nodeJ)/(1-dis);
                fixY(nodeI,nodeJ)=0;
                fixY(nodeJ,nodeI)=0;
            end
            % 线路三相短路故障处理结束
        end
        % 故障点修正完成
        % 节点导纳矩阵修正完成

        % 输出节点导纳矩阵
        %         fprintf(myf, '--------------Node Admittance Matrix----------\n');
        %         [count_i, count_j] = size(fixY);
        %         for i=1:count_i
        %             for k=1:count_j
        %                 fprintf(myf, ' %9.6f+j%10.6f ', real(fixY(i,k)),imag(fixY(i,k)));
        %             end
        %             fprintf(myf, '\n');
        %         end

        % 求取故障时的节点注入电流
        Iinj = zeros(nb,1);
        for k=1:nb
            busType=bus(k,7);
            if(calcSettings(1)==1)
                v = 1;
            else
                v = bus(k,2);
            end
            if (busType==1)
                % 计算发电机节点注入电流，近似认为内电势为1.0
                Iinj(bus(k,1),1) =  Iinj(bus(k,1),1) + v/j/bus(k,8);
            end
        end
        % 求取故障时的节点注入电流完成

        % 根据 U=Y ^ (-1) * I 求取节点电压；
        U = fixY \ Iinj;


        % 求取各支路及故障点的短路电流
        % 获取故障节点编号
        faultI=fault(1,1);
        faultJ=fault(1,2);
        dis=fault(1,3);
        faultNode = 0;
        if(faultI==0)
            faultNode = faultJ;
        end
        if(faultJ==0)
            faultNode = faultI;
        end
        if(dis==1)&&(faultI*faultJ~=0)
            faultNode = faultJ;
        end
        if(dis==0)&&(faultI*faultJ~=0)
            faultNode = faultI;
        end

        % 计算非故障支路的短路电流
        for k=1:nl
            I=line(k,1);
            J=line(k,2);
            Ui = U(I);
            if(J==0)
                Uj = 0;
            else
                Uj = U(J);
            end

            Ib(k,1) = I;
            Ib(k,2) = J;

            if J == 0
                Ym = line(k,5) + j * line(k,6);
                Ib(k,3) = Ui * Ym;
                Ib(k,4) = abs(Ib(k,3));
            else
                Zt=line(k,3)+j*line(k,4);
                Yt=1/Zt;
                Ym=line(k,5)+j*line(k,6);
                Ib(k,3) = Ui * (Ym + Yt) - Uj * Yt;
                Ib(k,4) = abs(Ib(k,3));
            end
        end

        % 节点三相短路，faultNode 不为 0
        if(faultNode~=0)
            % 前两列都显示故障节点编号
            Ib(nl+1,1) = faultNode;
            Ib(nl+1,2) = faultNode;
            Ifault = 0;
            for k=1:nl
                I=line(k,1);
                J=line(k,2);
                if(I*J==0)
                    continue;
                end
                Ib(k,1) = I;
                Ib(k,2) = J;
                % Ib(k,3)为已求得的非故障相关支路短路电流
                if (I==faultNode)
                    Ifault = Ifault - Ib(k,3);
                elseif (J==faultNode)
                    Ifault = Ifault + Ib(k,3);
                end
            end
            Ib(nl+1,3) = Ifault;
            Ib(nl+1,4) = abs(Ib(nl+1,3));
        end

        % 支路三相短路
        if(dis~=0)&&(dis~=1)&&(faultI*faultJ~=0)
            % 前两列分别表示故障支路的首末端节点编号
            for i = 1:nb
                if faultI == nodenum(i,2)
                    real_I = nodenum(i,1);
                end
                if faultJ == nodenum(i,2)
                    real_J = nodenum(i,1);
                end
            end
            Ib(nl+1,1) = real_I;
            Ib(nl+1,2) = real_J;
            Ib(nl+1,3) = U(real_I)*Y(real_I,real_J)/dis + U(real_J)*Y(real_I,real_J)/(1-dis);
            Ib(nl+1,4) = abs(Ib(nl+1,3));
        end

        [bus_new_result, line_result] = calculate_short_result (bus, Y, U, Ib, fault);


        %节点I   短路电流有效值
        [count_i, count_j] = size(line_result);
        if max_short_current < line_result(count_i,4)
            max_short_current = line_result(count_i,4);
            max_index_bus = index_bus;
            max_bus_new_result = bus_new_result;
            max_line_result = line_result;
        end
        fprintf(myf, '%d\t%9.6f\n', index_bus, line_result(count_i,4));
        fprintf(myf_graph, '%d\t%9.6f\n', index_bus, line_result(count_i,4));
    end
    
    bus_new_result = max_bus_new_result;
    line_result =max_line_result;

    [count_i, count_j] = size(line_result);
    fprintf(myf, 'Max short current Node and Value: \n');
    fprintf(myf, '%d\t%9.6f\n', max_index_bus, max_line_result(count_i,4));
    fprintf(myf_graph, 'Max short current Node and Value: \n');
    fprintf(myf_graph, '%d\t%9.6f\n', max_index_bus, max_line_result(count_i,4));
    
    [count_i, count_j] = size(bus_new_result);
    %节点电压结果
    fprintf(myf, 'Node Voltage calculation results: \n');
    fprintf(myf_graph, 'Node Voltage calculation results: \n');
    for i=1:count_i
        fprintf(myf, '%d\t', bus_new_result(i,1));
        fprintf(myf, '%9.6f+j%10.6f\t', real(bus_new_result(i,2)), imag(bus_new_result(i,2)));
        fprintf(myf, '%d\t', bus_new_result(i,3));
        fprintf(myf, '%d\t', bus_new_result(i,4));
        fprintf(myf, '\n');
        fprintf(myf_graph, '%d\t', bus_new_result(i,1));
        fprintf(myf_graph, '%9.6f+j%10.6f\t', real(bus_new_result(i,2)), imag(bus_new_result(i,2)));
        fprintf(myf_graph, '%d\t', bus_new_result(i,3));
        fprintf(myf_graph, '%d\t', bus_new_result(i,4));
        fprintf(myf_graph, '\n');
    end

    %各支路短路电流Ib
    fprintf(myf, 'Short circuit current of each branch Ib = \n');
    fprintf(myf_graph, 'Short circuit current of each branch Ib = \n');
    %节点I    节点J          短路电流相量    短路电流有效值
    fprintf(myf, 'Node I    Node J    Short circuit current phasor    Short circuit current effective value \n');
    fprintf(myf_graph, 'Node I    Node J    Short circuit current phasor    Short circuit current effective value \n');
    [count_i, count_j] = size(line_result);
    for i=1:count_i
        fprintf(myf, '%d\t', line_result(i,1));
        fprintf(myf, '%d\t', line_result(i,2));
        fprintf(myf, '%9.6f+j%10.6f\t', real(line_result(i,3)), imag(line_result(i,3)));
        fprintf(myf, '%9.6f\t', line_result(i,4));
        fprintf(myf, '\n');

        fprintf(myf_graph, '%d\t', line_result(i,1));
        fprintf(myf_graph, '%d\t', line_result(i,2));
        fprintf(myf_graph, '%9.6f+j%10.6f\t', real(line_result(i,3)), imag(line_result(i,3)));
        fprintf(myf_graph, '%9.6f\t', line_result(i,4));
        fprintf(myf_graph, '\n');
    end

    t=toc;
    fprintf(myf, 'Time for method: \n');
    fprintf(myf, num2str(t));
    fprintf(myf_graph, 'Time for method: \n');
    fprintf(myf_graph, num2str(t));
    fclose(myf);
    fclose(myf_graph);

end