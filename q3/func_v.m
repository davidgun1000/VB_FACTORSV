function [out]=func_v(A)
   %[m,n]=size(A); 
   %out=EliminationM(m)*vec(A);
    out=A(itril(size(A)));
end