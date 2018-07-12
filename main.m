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

%节点重新编号开始
[bus, line] = rearrange(bus, line);
%节点重新编号结束

Y= generateY(bus,line);
file_output('.\output.dat', '--------------节点导纳矩阵----------\n', Y);

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
    [deltaP,deltaQ] = dPQ(Y,bus);
    J=form_jac(Y,bus);
    
    deltaPQ=[deltaP;deltaQ];
    deltaUtheta = J^(-1)*deltaPQ;
    
    deltatheta = deltaUtheta(1:nPoint-nSW,:);
    deltaU = deltaUtheta(nPoint:nPoint+nPQ-nSW,:);
    
    theta(1:nPoint-nSW,1)=theta(1:nPoint-nSW,1) - deltatheta;
    U(1:nPQ,1)=U(1:nPQ,1) - deltaU;
    if (max(abs(deltatheta))<eps1)&&(max(abs(deltaU))<eps2)
        break;
    end
end
