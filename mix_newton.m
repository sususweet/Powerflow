clear;

max1=100;
eps1=1.0e-10;
eps2=1.0e-10;

[dfile,pathname]=uigetfile('*.m','Select Data File');
if pathname == 0
    error(' you must select a valid data file')
else
    lfile =length(dfile);
    % strip off .m
    eval(dfile(1:lfile-2));
end

global nodenum;
global nSW;
global nPV;
global nPQ;
global U;
global theta;
global P;
global Q;

%节点重新编号开始
[bus, line] = rearrange(bus, line);
%节点重新编号结束
filename = ['.\output_mix_newton_', dfile, '.dat'];
myf=fopen(filename,'w');

Y= generateY(bus,line);

fprintf(myf, '--------------节点导纳矩阵----------\n');
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
% 高斯-塞德尔法初值计算
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

    fprintf(myf, '\n');
    fprintf(myf, '-------------第%d次迭代的节点相角和电压的偏差dx----------\n', count);
    for i=1: length(deltatheta)
        fprintf(myf, 'dang%d   %13.6e\n', i, deltatheta(i,1));
    end
    for i=1: length(deltaU)
        fprintf(myf, 'dU%d     %13.6e\n', i, deltaU(i,1));
    end
    
    U = U_new;
    theta = theta_new;
    
    fprintf(myf, '\n');
    fprintf(myf, '-------------第%d次迭代的节点相角delta(弧度为单位）和电压U----------\n', count);
    for i=1: length(deltatheta)
        fprintf(myf, 'ang%d   %10.6f\n', i, theta(i,1));
    end
    for i=1: length(deltaU)
        fprintf(myf, 'U%d   %12.6f\n', i, U(i,1));
    end
% 高斯-塞德尔法初值计算结束

% 牛顿-拉夫逊法迭代开始
for count=1:max1
    fprintf(myf, '\n');
    fprintf(myf, '--------------第%d次迭代的结果----------\n', count);
    
    [deltaP,deltaQ] = dPQ(Y,bus);
    J=form_jac(Y,bus);
    
    fprintf(myf, '-------------第%d次迭代的雅可比矩阵J----------\n', count);
    [count_i, count_j] = size(J);
    for i=1:count_i
        for k=1:count_j
            fprintf(myf, '%11.6f  ', J(i,k));
        end
        fprintf(myf, '\n');
    end
    
    fprintf(myf, '\n');
    fprintf(myf, '-------------第%d次迭代的功率偏差dP和dQ----------\n', count);
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
    
    fprintf(myf, '\n');
    fprintf(myf, '-------------第%d次迭代的节点相角和电压的偏差dx----------\n', count);
    for i=1: length(deltatheta)
        fprintf(myf, 'dang%d   %13.6e\n', i, deltatheta(i,1));
    end
    for i=1: length(deltaU)
        fprintf(myf, 'dU%d/U%d     %13.6e\n', i, i, deltaU(i,1));
    end
    
    theta(1:nPoint-nSW,1)=theta(1:nPoint-nSW,1) - deltatheta;
    U(1:nPQ,1)=U(1:nPQ,1) - deltaU .* U(1:nPQ,1) ;
    
    fprintf(myf, '\n');
    fprintf(myf, '-------------第%d次迭代的节点相角delta(弧度为单位）和电压U----------\n', count);
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
fprintf(myf, '---------------混合牛顿－拉夫逊法潮流计算结果----------\n');
fprintf(myf, ' 节点计算结果：\n');
fprintf(myf, '节点   节点电压    节点相角（角度）    节点注入功率\n');

[count_i, count_j] = size(bus_new_result);
for i=1:count_i
    fprintf(myf, ' %d  ', bus_new_result(i,1));
    for k=2:3
        fprintf(myf, '%10.6f  ', bus_new_result(i,k));
        fprintf(myf, '   ');
    end
    fprintf(myf, '%9.6f+j%10.6f ', bus_new_result(i,4), bus_new_result(i,5));
    fprintf(myf, '\n');
end

fprintf(myf, '\n');
fprintf(myf, ' 线路计算结果：\n');
fprintf(myf, '节点I    节点J          线路功率S(I,J)         线路功率S(J,I)         线路损耗dS(I,J)\n');
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
end
fclose(myf);
