function status = file_output(filename, title, matrix)
myf=fopen(filename,'w');
fprintf(myf, title);
[count_i, count_j] = size(matrix);

for i=1:count_i
    for j=1:count_j
        fprintf(myf, ' %8.6f+j %8.6f', real(matrix(i,j)),imag(matrix(i,j)));
    end
    fprintf(myf, '\n');
end

fclose(myf);
status = 1;