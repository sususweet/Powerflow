Sb=100; 
Ub2=121;
Ub1=10.5;
Ub=110;
Zb=Ub^2/Sb;
Yb=1/Zb;

bus = [ 1 1.00    0.00  -0.35  -0.18      1;
	    2 1.00    0.00  -0.30  -0.16      1;
	    3 1.00    0.00  -0.40  -0.20      1;
	    4 1.00    0.00  -0.25  -0.12      1;
	    5 1.00    0.00  -0.20  -0.10      1;
	    6 1.00    0.00  -0.15  -0.06      1;
	    7 1.05    0.00   0.40   0.00      2;
	    8 1.05    0.00   0.45   0.00      2;
	    9 1.05    0.00   0.00   0.00      3];

% line data format
% line: from bus, to bus, resistance(pu), reactance(pu),
%       line charging(pu), tap ratio

line = [ 2 5  0.131   0.400      0   2.98e-6                  0 25;
         1 4  0.131   0.400      0   2.98e-6                  0 50;
         1 6  0.105   0.393      0   3.04e-6                  0 45;
         3 5  0.131   0.400      0   2.98e-6                  0 30;
         3 6  0.105   0.393      0   3.04e-6                  0 35;
         2 4  0.131   0.400      0   2.98e-6                  0 20;
         4 7  250     10.5      55     1.5                  1.0 50;
         6 9  410     10.5      98     1.5                  1.0 80;
         5 8  250     10.5      55     1.5                  1.0 50];

[nnl,mml]=size(line);
for i=1:nnl,
    if line(i,7)==0
        line(i,3)=line(i,3)*line(i,8)/Zb;
        line(i,4)=line(i,4)*line(i,8)/Zb;
        line(i,5)=line(i,5)*line(i,8)/Yb;
        line(i,6)=line(i,6)*line(i,8)/2/Yb;
    else
        Pk=line(i,3);
        Uk=line(i,4);
        P0=line(i,5);
        I0=line(i,6);
        Sn=line(i,8);
        line(i,3)=Pk/1000/Sn/Sn*Ub2*Ub2/Zb;
        line(i,4)=Uk/100/Sn*Ub2*Ub2/Zb;
        line(i,5)=P0/1000/Ub2/Ub2/Yb;
        line(i,6)=-I0/100*Sn/Ub2/Ub2/Yb;
        line(i,7)=(10.5/121)/(10/110);
    end
end

line=line(:,1:7);