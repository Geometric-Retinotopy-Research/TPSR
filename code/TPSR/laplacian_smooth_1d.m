function z = laplacian_smooth_1d(y,w,lambda)
m = length(y);
E = eye(m);
D = diff(E);
z =  (diag(w) + lambda * (D' * D))\ (diag(w)*y);
end