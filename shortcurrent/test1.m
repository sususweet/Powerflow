% 三相短路电流计算测试系统
% 平启动时无需考虑bus矩阵末两列
% data for test case

% 节点参数
% (bus#)( volt )( ang )(  p )(  q )(bus type)(是否机组1是0否)(Xd'')(calc voltage)(calc ang)

bus=[
   1   1.00  0.0   -1.60  -0.80   1  0   0.000  0.862150  -4.778511;
   2   1.00  0.0   -2.00  -1.00   1  0   0.000  1.077916  17.853530;
   3   1.00  0.0   -3.70  -1.30   1  0   0.000  1.036411  -4.281930;
   4   1.05  0.0    5.00   0.00   2  1   0.020  1.050000  21.843319;
   5   1.05  0.0    0.00   0.00   3  1   0.020  1.050000  0.000000];

% 支路参数
% b#1 b#2  ( R )( X )( G )( B/2 ) ( tap )
line = [
   1   2 0.04 0.25  0.0 0.25    0;
   1   3 0.10 0.35  0.0 0.0     0;
   2   3 0.08 0.30  0.0 0.25    0;
   5   3 0.00 0.03  0.0 0.0  1.05;
   4   2 0.00 0.015 0.0 0.0  1.05];

% 故障参数
% 故障支路首端节点 故障支路末端短路 （距首端距离0-1）
fault = [
    1   3   0.5];
   

% 计算参数
% 是否平启动（1为是，0为否）
calcSettings = [
    1];

disp('测试系统：双机五节点测试系统');
disp('故障类型：三相短路故障');
disp('故障点：支路1-3正中');
disp('短路电阻：0');
disp('计算类型：平启动');
    
