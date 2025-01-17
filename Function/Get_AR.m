function  AR  =  Get_AR( X, ws )

b2        =  size(X, 1);
b         =  sqrt(b2);
c         =  ceil(b/2);
r         =  floor(ws/2);

M         =  zeros(b, b);
M(c-r:c+r,c-r:c+r)  =  1;
M         =  M(:);
ind       =  find( M );

X         =  X(ind, :);
b2        =  size(X, 1);

cen       =  ceil( b2/2 );
Y         =  X(cen,:)';
X(cen,:)  =  [];
X         =  X';

A         =   X'*X;
b         =   X'*Y;
wei       =  lsqminnorm(A,b);
err1      =  mean( (Y-X*wei).^2 );


if  ws>3
    
    wei      =  l1_ls(X./255, Y./255, 0.01, 0.001);
    err2     =  mean( (Y-X*wei).^2 );
        
end


AR       =  wei;
