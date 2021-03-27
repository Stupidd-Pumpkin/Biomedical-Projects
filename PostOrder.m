function matrix= PostOrder(matrix,index,x)
if matrix(index,3)
    matrix=PostOrder(matrix,matrix(index,3),2*x+1);
    matrix=PostOrder(matrix,matrix(index,4),2*x);
end
matrix(index,5)=x;
end