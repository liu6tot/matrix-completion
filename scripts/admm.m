function [X, U, Y, Z, A] = admm(U_i, Y_i, Z_i, M, Omega, l1_norm)
   pho = 10^8;
   A = 0;
   tol = 10^-4;
   max_iter = 100;
   U = U_i;
   Y = Y_i;
   Z = Z_i;
   gamma = 1;
   t = 1;
   output_size = size(M, 1);
   
   while(true)
       % Update X 
       UY = U*Y;
       X = (Z*(UY'))/(UY*(UY'));
       
       % Update Y
       XU = X*U;
       Y = ((XU')*(XU))\((XU')*Z);
       
       % Update U
       U = pinv(X'*X)*(X'*Z*Y')*pinv(Y*Y');
       
       % Update Z
       if l1_norm == 1
          Z = X*U*Y + Omega.*(M - X*U*Y) - A/pho;

          D = DCT_Matrix();
          Dt = D';
          y = Z(:);

          m = size(y, 1);
          n = m;

          lambda  = 0.01; % regularization parameter
          rel_tol = 0.01; % relative target duality gap
          quiet = 1;

          %run the l1-regularized least squares solver
          [x, ~]=l1_ls(D,Dt,m,n,y,lambda,rel_tol, quiet);

          Z_opt = D*x;
          Z = reshape(Z_opt, [output_size, output_size]);
       else
          Z = X*U*Y + Omega.*(M - X*U*Y) - A/pho;
       end
       
       % Update A
       A = A + gamma*pho*Omega.*(Z - M);
       
       % Update iteration count.
       t = t + 1;
       
       disp(t);
       
       if t == max_iter
           break
       end
    end
end