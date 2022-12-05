function e_TSVD=tsvd(U,D,V,lambda,max_s)
D_diag=diag(D);
B=abs(D_diag-lambda*max_s);
[x,index]=sort(B);
cc=index(1);
D_diag_i=1./D_diag;
D_diag_i_t=D_diag_i(1:cc);
D_i_t=diag(D_diag_i_t);

U_t=U';
U_t=U_t(1:cc,:);
V=V(:,1:cc);
e_TSVD=V*D_i_t*U_t;
fprintf('  --> Truncated %.1f percent of singular values (cc: %d)\n',(numel(D_diag)-cc)/numel(D_diag)*100,cc);
end