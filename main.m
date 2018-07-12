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

global nSW;
global nPV;
global nPQ;
global U;
global theta;
global P;
global Q;

%�ڵ����±�ſ�ʼ
[bus, line] = rearrange(bus, line);
%�ڵ����±�Ž���

myf=fopen('.\output.dat','w');

Y= generateY(bus,line);

fprintf(myf, '--------------�ڵ㵼�ɾ���----------\n');
[count_i, count_j] = size(Y);
for i=1:count_i
    for j=1:count_j
        fprintf(myf, ' %9.6f+j%10.6f ', real(Y(i,j)),imag(Y(i,j)));
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
    fprintf(myf, '--------------��%d�ε����Ľ��----------\n', count);
    
    [deltaP,deltaQ] = dPQ(Y,bus);
    J=form_jac(Y,bus);
    
    fprintf(myf, '-------------��%d�ε������ſɱȾ���J----------\n', count);
    [count_i, count_j] = size(J);
    for i=1:count_i
        for j=1:count_j
            fprintf(myf, '%11.6f  ', J(i,j));
        end
        fprintf(myf, '\n');
    end
    
    fprintf(myf, '\n');
    fprintf(myf, '-------------��%d�ε����Ĺ���ƫ��dP��dQ----------\n', count);
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
    fprintf(myf, '-------------��%d�ε����Ľڵ���Ǻ͵�ѹ��ƫ��dx----------\n', count);
    for i=1: length(deltatheta)
        fprintf(myf, 'dang%d   %13.6e\n', i, deltatheta(i,1));
    end
	for i=1: length(deltaU)
        fprintf(myf, 'dU%d     %13.6e\n', i, deltaU(i,1));
    end
    
    theta(1:nPoint-nSW,1)=theta(1:nPoint-nSW,1) - deltatheta;
    U(1:nPQ,1)=U(1:nPQ,1) - deltaU;
    
    fprintf(myf, '\n');
    fprintf(myf, '-------------��%d�ε����Ľڵ����delta(����Ϊ��λ���͵�ѹU----------\n', count);
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
fclose(myf);

