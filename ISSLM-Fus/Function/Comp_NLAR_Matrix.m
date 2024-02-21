function  W1   =  Comp_NLAR_Matrix( Y, sz )
Y       =   Y*255;
ch      =   size(Y,1);
h       =   sz(1);
w       =   sz(2);
im      =   zeros(h, w, ch);
for i   =   1 : ch
    im(:,:,i)     =  reshape( Y(i,:), h, w );
end
S          =   18;
f          =   3;
f2         =   f^2; 
t          =   floor(f/2);
nv         =   10;
e_im       =   padarray( im, [t t], 'symmetric' );
[h w]      =   size( im(:,:,1) );
nt         =   (nv)*h*w;
R          =   zeros(nt,1);
C          =   zeros(nt,1);
V          =   zeros(nt,1);
L          =   h*w;
X          =   zeros(f*f*ch, L, 'single');
k          =   0;
lam_w      =  48000;
lamada     =  0.72;

for i  = 1:f
    for j  = 1:f
        
        k        =   k+1;
        blk      =   e_im(i:end-f+i,j:end-f+j, 1);
        X(k,:)   =   blk(:)'; 
        
        if ch>1
            blk  =  e_im(i:end-f+i,j:end-f+j, 2);
            X(k+f2,:) =  blk(:)';
            if ch>2
                blk  =  e_im(i:end-f+i,j:end-f+j, 3);
                X(k+f2*2,:) =  blk(:)';
            end            
        end 
     
    end
end
X           =   X'; 
X2          =   sum(X.^2, 2);
f2          =   f^2;
I           =   reshape((1:L), h, w);
f3          =   f2*ch;
cnt         =  1;
for  row  =  1 : h
    for  col  =  1 : w
        
        off_cen  =  (col-1)*h + row;        
        
        rmin    =   max( row-S, 1 );
        rmax    =   min( row+S, h );
        cmin    =   max( col-S, 1 );
        cmax    =   min( col+S, w );
         
        idx     =   I(rmin:rmax, cmin:cmax);
        idx     =   idx(:);
        B       =   X(idx, :);      
        B2      =   X2(idx, :);
        v       =   X(off_cen, :);
        v2      =   X2(off_cen, :);
        c2      =   B*v';
        
        dis     =   (B2 + v2 - 2*c2)/f3;
        [val,ind]     =   sort(dis);   

        b       =   B( ind(2:nv + 1),: )*v';
        wei     =   cgsolve(B( ind(2:nv + 1),: )*B( ind(2:nv + 1),: )' + lam_w*eye(nv), b); 
        
        R(cnt:cnt+nv)     =   off_cen; 
        C(cnt:cnt+nv)     =   [idx( ind(2:nv+1) );  off_cen];
        V(cnt:cnt+nv)     =   [lamada*(wei./(sum(wei)+eps)); (1-lamada)];  
       
        cnt                 =   cnt + nv + 1;        
    end
end
R     =   R(1:cnt-1);
C     =   C(1:cnt-1);
V     =   V(1:cnt-1);
W1    =   sparse(R, C, V, h*w, h*w);
W1    =   W1';  
