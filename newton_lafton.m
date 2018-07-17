function Y = newton_lafton(dfile)
    max1=100;
    eps1=1.0e-10;
    eps2=1.0e-10;

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

        [deltaP,deltaQ] = dPQ(Y,bus);
        J=form_jac(Y,bus);

        %第x次迭代的雅可比矩阵J
        fprintf(myf, '-------------The Jacobi matrix J of the %d iteration----------\n', count);
        [count_i, count_j] = size(J);
        for i=1:count_i
            for k=1:count_j
                fprintf(myf, '%11.6f  ', J(i,k));
            end
            fprintf(myf, '\n');
        end

        %第x次迭代的功率偏差dP和dQ
        fprintf(myf, '\n');
        fprintf(myf, '-------------The power deviation dP and dQ of the %d iteration----------\n', count);
        for i=1: length(deltaP)
            fprintf(myf, 'dP%d   %13.6e\n', i, deltaP(i,1));
        end
        for i=1: length(deltaQ)
            fprintf(myf, 'dQ%d   %13.6e\n', i, deltaQ(i,1));
        end

        deltaPQ=[deltaP;deltaQ];
        deltaUtheta = J^(-1)*deltaPQ;

        deltatheta = deltaUtheta(1:nPoint-nSW,:);
        deltaU = deltaUtheta(nPoint:nPoint+nPQ-nSW,:);

        %第x次迭代的节点相角和电压的偏差dx
        fprintf(myf, '\n');
        fprintf(myf, '-------------The node phase angle and voltage deviation dX of the %d iteration----------\n', count);
        for i=1: length(deltatheta)
            fprintf(myf, 'dang%d   %13.6e\n', i, deltatheta(i,1));
        end
        for i=1: length(deltaU)
            fprintf(myf, 'dU%d/U%d     %13.6e\n', i, i, deltaU(i,1));
        end
        fprintf(myf_iteration, '%d\t%13.6e\t%13.6e\n', count, max(abs(deltatheta)), max(abs(deltaU .* U(1:nPQ,1))));

        theta(1:nPoint-nSW,1)=theta(1:nPoint-nSW,1) - deltatheta;
        U(1:nPQ,1)=U(1:nPQ,1) - deltaU .* U(1:nPQ,1) ;
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
    %牛顿－拉夫逊法潮流计算结果
    fprintf(myf, '---------------The results of the Newton LV flow calculation----------\n');
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