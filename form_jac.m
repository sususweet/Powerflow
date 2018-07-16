function J=form_jac(Y,bus)
    global nPQ;
    global nPV;
    global nSW;
    global U;
    global theta;

    nPoint = length(Y);
    G=real(Y);
    B=imag(Y);

    %开始计算雅可比矩阵
    for row=1:nPoint-nSW
        for col = 1:nPoint-nSW
                if row == col
                    P_cache = 0;
                    for loop = 1:nPoint
                        P_cache = P_cache + U(row)*U(loop)*(G(row,loop)*cos(theta(row)-theta(loop))+B(row,loop)*sin(theta(row)-theta(loop)));
                    end
                    Q_cache = 0;
                    for loop = 1:nPoint
                        Q_cache = Q_cache + U(row)*U(loop)*(G(row,loop)*sin(theta(row)-theta(loop))-B(row,loop)*cos(theta(row)-theta(loop)));
                    end

                    H(row,row) = U(row)^2*B(row,row)+Q_cache;
                    if (row <= nPQ)
                        N(row,row) = -U(row)^2*G(row,row)-P_cache;
                        J(row,row)=N(row,row)+2*U(row)^2*G(row,row);
                        L(row,row) = -H(row,row) +2*U(row)^2*B(row,row);
                    end
                else
                    H(row,col) = -U(row) * U(col) *(G(row,col) * sin(theta(row)-theta(col)) - B(row,col) * cos(theta(row)-theta(col)));
                    if (col <= nPQ)
                        N(row,col) = -U(row) * U(col) *(G(row,col) * cos(theta(row)-theta(col)) + B(row,col) * sin(theta(row)-theta(col)));
                    end
                    if (row <= nPQ)
                        J(row,col) = -(-U(row) * U(col) *(G(row,col) * cos(theta(row)-theta(col)) + B(row,col) * sin(theta(row)-theta(col))));
                    end
                    if (row <= nPQ && col<= nPQ)
                        L(row,col) = -U(row) * U(col) *(G(row,col) * sin(theta(row)-theta(col)) - B(row,col) * cos(theta(row)-theta(col)));
                    end
                end
        end
    end
    %雅可比矩阵计算完成
    J=[H,N;J,L];
end