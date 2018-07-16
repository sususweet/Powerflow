% data for New England Test case
disp('New England data')
%(bus#)( volt )( ang  )(  p )(  q )(bus type)(是否机组1是0否)(Xd'')
bus=[
   1   1.0000   0.00   0.0000   0.0000   1  0  0.000;
   2   1.0000   0.00   0.0000   0.0000   1  0  0.000;
   3   1.0000   0.00  -3.2200  -0.0240   1  0  0.000;
   4   1.0000   0.00  -5.0000  -1.8400   1  0  0.000;
   5   1.0000   0.00   0.0000   0.0000   1  0  0.000;
   6   1.0000   0.00   0.0000   0.0000   1  0  0.000;
   7   1.0000   0.00  -2.3380  -0.8400   1  0  0.000;
   8   1.0000   0.00  -5.2200  -1.7600   1  0  0.000;
   9   1.0000   0.00   0.0000   0.0000   1  0  0.000;
  10   1.0000   0.00   0.0000   0.0000   1  0  0.000;
  11   1.0000   0.00   0.0000   0.0000   1  0  0.000;
  12   1.0000   0.00  -0.0850  -0.8800   1  0  0.000;
  13   1.0000   0.00   0.0000   0.0000   1  0  0.000;
  14   1.0000   0.00   0.0000   0.0000   1  0  0.000;
  15   1.0000   0.00  -3.2000  -1.5300   1  0  0.000;
  16   1.0000   0.00  -3.2940  -0.3230   1  0  0.000;
  17   1.0000   0.00   0.0000   0.0000   1  0  0.000;
  18   1.0000   0.00  -1.5800  -0.3000   1  0  0.000;
  19   1.0000   0.00   0.0000   0.0000   1  0  0.000;
  20   1.0000   0.00  -6.8000  -1.0300   1  0  0.000;
  21   1.0000   0.00  -2.7400  -1.1500   1  0  0.000;
  22   1.0000   0.00   0.0000   0.0000   1  0  0.000;
  23   1.0000   0.00  -2.4750  -0.8460   1  0  0.000;
  24   1.0000   0.00  -3.0800   0.9220   1  0  0.000;
  25   1.0000   0.00  -2.2400  -0.4720   1  0  0.000;
  26   1.0000   0.00  -1.3900  -0.1700   1  0  0.000;
  27   1.0000   0.00  -2.8100  -0.7550   1  0  0.000;
  28   1.0000   0.00  -2.0600  -0.2760   1  0  0.000;
  29   1.0000   0.00  -2.8350  -0.2690   1  0  0.000;
  30   1.0475   0.00   2.5000   0.0000   2  1  0.0279;
  31   1.0000   0.00   0.0000   0.0000   3  1  0.0627;
  32   1.0000   0.00   6.5000   1.7590   1  1  0.0478;
  33   1.0000   0.00   6.3200   1.0335   1  1  0.0392;
  34   1.0123   0.00   5.0800   0.0000   2  1  0.1188;
  35   1.0493   0.00   6.5000   0.0000   2  1  0.045;
  36   1.0000   0.00   5.6000   0.9688   1  1  0.0441;
  37   1.0278   0.00   5.4000   0.0000   2  1  0.0513;
  38   1.0265   0.00   8.3000   0.0000   2  1  0.0513;
  39   1.0300   0.00  -1.0400   0.0000   2  1  0.0128];

%b#1 b#2( res  )(  rea  )(charg )(tap )(phase)
line = [
   1   2 0.00350  0.04110 0 0.34935 0;
   1  39 0.00100  0.02500 0 0.3750  0;
   2   3 0.00130  0.01510 0 0.12860 0;
   2  25 0.00700  0.00860 0 0.0730  0;
   3   4 0.00130  0.02130 0 0.11070 0;
   3  18 0.00110  0.01330 0 0.10690 0;
   4   5 0.00080  0.01280 0 0.06710 0;
   4  14 0.00080  0.01290 0 0.06910 0;
   6   5 0.00020  0.00260 0 0.02170 0;
   5   8 0.00080  0.01120 0 0.07380 0;
   6   7 0.00060  0.00920 0 0.05650 0;
   6  11 0.00070  0.00820 0 0.06945 0;
   7   8 0.00040  0.00460 0 0.03900 0;
   8   9 0.00230  0.03630 0 0.19020 0;
   9  39 0.00100  0.02500 0 0.60000 0;
  10  11 0.00040  0.00430 0 0.03645 0;
  10  13 0.00040  0.00430 0 0.03645 0;
  13  14 0.00090  0.01010 0 0.08615 0;
  14  15 0.00180  0.02170 0 0.18300 0;
  15  16 0.00090  0.00940 0 0.08550 0;
  16  17 0.00070  0.00890 0 0.06710 0;
  16  19 0.00160  0.01950 0 0.15200 0;
  16  21 0.00080  0.01350 0 0.12740 0;
  16  24 0.00030  0.00590 0 0.03400 0;
  17  18 0.00070  0.00820 0 0.06595 0;
  17  27 0.00130  0.01730 0 0.16080 0;
  21  22 0.00080  0.01400 0 0.12825 0;
  22  23 0.00060  0.00960 0 0.09230 0;
  23  24 0.00220  0.03500 0 0.18050 0;
  25  26 0.00320  0.03230 0 0.25650 0;
  26  27 0.00140  0.01470 0 0.11980 0;
  26  28 0.00430  0.04740 0 0.39010 0;
  26  29 0.00570  0.06250 0 0.51450 0;
  28  29 0.00140  0.01510 0 0.12450 0;
   4   0 0.00000  0.00000 0 1.00000 0;
   5   0 0.00000  0.00000 0 2.00000 0;
  11  12 0.00160  0.04350 0 0.00000 1.0060;
  13  12 0.00160  0.04350 0 0.00000 1.0060;
  30   2 0.00000  0.01810 0 0.00000 1.0250;
  31   6 0.00000  0.02500 0 0.00000 1.0700;
  32  10 0.00000  0.02000 0 0.00000 1.0700;
  34  20 0.00090  0.01800 0 0.00000 1.0090;
  33  19 0.00070  0.01420 0 0.00000 1.0700;
  35  22 0.00000  0.01430 0 0.00000 1.0250;
  36  23 0.00050  0.02720 0 0.00000 1.0000;
  37  25 0.00060  0.02320 0 0.00000 1.0250;
  38  29 0.00080  0.01560 0 0.00000 1.0250;
  20  19 0.00070  0.01380 0 0.00000 1.0600  
];


% 故障参数
% 故障支路首端节点 故障支路末端短路 （距首端距离0-1）
fault = [
   22   0   0];

% 计算参数
% 是否平启动（1为是，0为否）
calcSettings = [
    1];

disp('测试系统：新英格兰-十机三十九节点测试系统');
disp('故障类型：三相短路故障');
disp('故障点：节点22');
disp('短路电阻：0');
disp('计算类型：平启动');