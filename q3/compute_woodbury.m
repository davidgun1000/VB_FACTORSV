function [result]=compute_woodbury(B,d)
num_factor=size(B,2);
D_square_inv=diag([1./(d.^2)]);
result=D_square_inv-D_square_inv*B*((eye(num_factor)+B'*D_square_inv*B)\(B'*D_square_inv));
end