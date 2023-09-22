function varargout=sqrth(A)
% compute decomposition A=F'*F and inv(F) 
% using singular value decomposition
%
% If you want to decompose an hermitian matrix
% as A=F*F' (covariance stuff)
% you can do F=sqrth(A)' ; and check with norm(A-F*F')

[Uexp,Sexp,~]=svd(A);
s=sqrt( diag(Sexp) );
F=diag(s)*Uexp';
Finv=Uexp*diag(ones(size(s))./s);

switch nargout
 case 1
  varargout{1}=F;
 case 2
  varargout{1}=F;
  varargout{2}=Finv;
 otherwise
end
