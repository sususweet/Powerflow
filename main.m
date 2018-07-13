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
filename = ['.\output_', dfile, '.dat'];
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
fprintf(myf, '\n');
fprintf(myf, '---------------牛顿－拉夫逊法潮流计算结果----------\n');
fprintf(myf, ' 节点计算结果：\n');
fprintf(myf, '节点   节点电压    节点相角（角度）    节点注入功率\n');

G=real(Y);
B=imag(Y);

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
