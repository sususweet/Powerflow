function Y = gassion_sader(dfile)
    max1=1000;
    eps1=1.0e-10;
    eps2=1.0e-10;

    lfile =length(dfile);
    % strip off .m
    eval(dfile(1:lfile-2));

    % [dfile,pathname]=uigetfile('*.m','Select Data File');
    % if pathname == 0
    %     error(' you must select a valid data file')
    % else
    %     lfile =length(dfile);
    %     % strip off .m
    %     eval(dfile(1:lfile-2));
    % end

    global nodenum;
    global nSW;
    global nPV;
    global nPQ;
    global U;
    global theta;
    global P;
    global Q;
    tic;

    %节点重新编号开始
    [bus, line] = rearrange(bus, line);
    %节点重新编号结束
    filename = ['./result/output_', dfile, '.dat'];
    myf=fopen(filename,'w','n','UTF-8');

    filename_graph = ['./result/output_graph_', dfile, '.dat'];
    myf_graph=fopen(filename_graph,'w','n','UTF-8');

    filename_iteration = ['./result/output_iteration_', dfile, '.dat'];
    myf_iteration=fopen(filename_iteration,'w','n','UTF-8');
    
    Y= generateY(bus,line);
    %输出节点导纳矩阵
    fprintf(myf, '--------------Node Admittance Matrix----------\n');
    [count_i, count_j] = size(Y);
    for i=1:count_i
        for k=1:count_j
            fprintf(myf, ' %9.6f+j%10.6f ', real(Y(i,k)),imag(Y(i,k)));
        end
        fprintf(myf, '\n');
    end

    nPoint = length(Y);

    U=zeros(nPoint,1);
    theta=zeros(nPoint,1);
    for i=1:nPoint
        U(i)=bus(i,2);
        theta(i)=bus(i,3);
        P(i,1) = bus(i,4);
        Q(i,1) = bus(i,5);
    end

    count = 0;
    fprintf(myf_iteration, '-------------The node phase angle and voltage deviation dX of the iteration----------\n');
    fprintf(myf_iteration, 'Iteration count      dang      dU\n');
    for count=1:max1
        fprintf(myf, '\n');
        %第x次迭代的结果
        fprintf(myf, '--------------The results of the %d iteration----------\n', count);

        U_new=zeros(nPoint,1);
        theta_new=zeros(nPoint,1);

        %采用异步迭代
        %首先处理PQ节点
        for i = 1:nPQ
            Ui = U(i,1) * exp(j*theta(i,1));
            sumYU = 0;
            for k = 1:nPoint
                if k < i
                    sumYU = sumYU + Y(i,k)*U_new(k,1) * exp(j*theta_new(k,1));
                elseif k > i
                    sumYU = sumYU + Y(i,k)*U(k,1) * exp(j*theta(k,1));
                end
            end
            Ui_new = ((P(i,1) - j*Q(i,1))/conj(Ui) - sumYU) / Y(i,i);

            U_new(i,1) = abs(Ui_new);
            theta_new(i,1) = angle(Ui_new);
        end
        %PQ节点处理完毕

        %PV节点处理开始
        for i = nPQ +1:nPQ+nPV
            Ui = U(i,1) * exp(j*theta(i,1));
            sumYU = 0;
            for k = 1:nPoint
                sumYU = sumYU + conj(Y(i,k))*conj(U(k,1) * exp(j*theta(k,1)));
            end
            Qp_new = imag(Ui*sumYU);

            Ui = U(i,1) * exp(j*theta(i,1));
            sumYU = 0;
            for k = 1:nPoint
                if k < i
                    sumYU = sumYU + Y(i,k)*U_new(k,1) * exp(j*theta_new(k,1));
                elseif k > i
                    sumYU = sumYU + Y(i,k)*U(k,1) * exp(j*theta(k,1));
                end
            end
            Ui_new = ((P(i,1) - j*Qp_new)/conj(Ui) - sumYU) / Y(i,i);

            %只对相角进行修正
            U_new(i,1) = U(i,1);
            theta_new(i,1) = angle(Ui_new);
        end
        %PV节点处理完毕

        %SW节点处理开始
        for i = nPQ+nPV +1:nPoint
            U_new(i,1) = U(i,1);
            theta_new(i,1) = theta(i,1);
        end
        %SW节点处理结束

        %计算误差值
        deltaU = U_new - U;
        deltatheta = theta_new - theta;

        %第x次迭代的节点相角和电压的偏差dx
        fprintf(myf, '\n');
        fprintf(myf, '-------------The node phase angle and voltage deviation dX of the %d iteration----------\n', count);
        for i=1: length(deltatheta)
            fprintf(myf, 'dang%d   %13.6e\n', i, deltatheta(i,1));
        end
        for i=1: length(deltaU)
            fprintf(myf, 'dU%d     %13.6e\n', i, deltaU(i,1));
        end
        fprintf(myf_iteration, '%d\t%13.6e\t%13.6e\n', count, max(abs(deltatheta)), max(abs(deltaU)));

		U = U_new;
        theta = theta_new;
        %第x次迭代的节点相角delta(弧度为单位）和电压U
        fprintf(myf, '\n');
        fprintf(myf, '-------------The node phase Delta (radian unit) and voltage U of the %d iteration----------\n', count);
        for i=1: length(deltatheta)
            fprintf(myf, 'ang%d   %10.6f\n', i, theta(i,1));
        end
        for i=1: length(deltaU)
            fprintf(myf, 'U%d   %12.6f\n', i, U(i,1));
        end

        if (max(abs(deltatheta))<eps1)&&(max(abs(deltaU))<eps2)
            break;
        end
    end

    [bus_new_result, line_result] = calculate_result (Y, bus, line);

    fprintf(myf, '\n');
    %高斯 - 塞德尔法潮流计算结果
    fprintf(myf, '---------------The result of the Gauss Seidel method calculation----------\n');
    %节点计算结果
    fprintf(myf, 'Node calculation results: \n');
    fprintf(myf_graph, 'Node calculation results: \n');
    %节点   节点电压    节点相角（角度）    节点注入功率    节点类型
    fprintf(myf, 'Node   Node voltage   Angle(degree)   Node injection power    Type\n');
    fprintf(myf_graph, 'Node   Node voltage   Angle(degree)   Node injection power    Type\n');

    [count_i, count_j] = size(bus_new_result);
    for i=1:count_i
        fprintf(myf, ' %d  ', bus_new_result(i,1));
        fprintf(myf_graph, '%d\t', bus_new_result(i,1));
        for k=2:3
            fprintf(myf, '%10.6f  ', bus_new_result(i,k));
            fprintf(myf, '   ');
            fprintf(myf_graph, '%10.6f\t', bus_new_result(i,k));
        end
        fprintf(myf, '%9.6f+j%10.6f\t', bus_new_result(i,4), bus_new_result(i,5));
        fprintf(myf, '%d', bus_new_result(i,6));
        fprintf(myf, '\n');
        fprintf(myf_graph, '%9.6f+j%10.6f\t', bus_new_result(i,4), bus_new_result(i,5));
        fprintf(myf_graph, '%d', bus_new_result(i,6));
        fprintf(myf_graph, '\n');
    end

    fprintf(myf, '\n');
    %线路计算结果
    fprintf(myf, 'Line calculation results: \n');
    fprintf(myf_graph, 'Line calculation results: \n');
    %节点I    节点J          线路功率S(I,J)         线路功率S(J,I)         线路损耗dS(I,J)
    fprintf(myf, 'Node I   Node J    Line power S(I, J)    Line power S(J, I)    Line loss dS(I, J)\n');
    fprintf(myf_graph, 'Node I   Node J    Line power S(I, J)    Line power S(J, I)    Line loss dS(I, J)\n');
    [count_i, count_j] = size(line_result);
    for i=1:count_i
        fprintf(myf, ' ');
        fprintf(myf, '%d', line_result(i,1));
        fprintf(myf, '         ');
        fprintf(myf, '%d', line_result(i,2));
        fprintf(myf, '         ');
        fprintf(myf, '%9.6f+j%10.6f', real(line_result(i,3)), imag(line_result(i,3)));
        fprintf(myf, '  ');
        fprintf(myf, '%9.6f+j%10.6f', real(line_result(i,4)), imag(line_result(i,4)));
        fprintf(myf, '  ');
        fprintf(myf, '%9.6f+j%10.6f', real(line_result(i,5)), imag(line_result(i,5)));
        fprintf(myf, ' \n');


        fprintf(myf_graph, '%d\t', line_result(i,1));
        fprintf(myf_graph, '%d\t', line_result(i,2));
        fprintf(myf_graph, '%9.6f+j%10.6f\t', real(line_result(i,3)), imag(line_result(i,3)));
        fprintf(myf_graph, '%9.6f+j%10.6f\t', real(line_result(i,4)), imag(line_result(i,4)));
        fprintf(myf_graph, '%9.6f+j%10.6f', real(line_result(i,5)), imag(line_result(i,5)));
        fprintf(myf_graph, '\n');
    end
    t=toc;
    fprintf(myf, 'Time for method: \n');
    fprintf(myf, num2str(t));
    fprintf(myf_graph, 'Time for method: \n');
    fprintf(myf_graph, num2str(t));
    fclose(myf);
    fclose(myf_graph);
    fclose(myf_iteration);
end