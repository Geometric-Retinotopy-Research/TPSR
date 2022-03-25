function [e_x,e_y]=grad(A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using Fast Marching Algorithm (FMA) one can obtain potential filed.  %
% In order to obtain vector field of desired direction one must find   %
% gradient of previously obtained potenial field. Since wall positions %
% in potential field have value of inf, it is not possible to use      %
% MATLAB function gradient. For this reason the following function was %
% was written and it allows to cope with the problem in effective way.
% A-input matrix from gradientmap.f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a b]=size(A);
e_x=zeros(a,b);
e_y=zeros(a,b);

for m=1:a
    for n=1:b
        
        %X-direction
        if (A(m,n)~=inf && n>1 && n<b)
           if (A(m,n-1))==inf 
               e_x(m,n)=A(m,n)-A(m,n+1);
           elseif (A(m,n+1)==inf)
               e_x(m,n)=A(m,n-1)-A(m,n);
           else
               e_x(m,n)=(A(m,n-1)-A(m,n+1))/2;
           end
        end
            
        if (A(m,n)==inf && n>1 && n<b)
           if (A(m,n-1)==inf || A(m,n+1)==inf)
               e_x(m,n)=0;
           elseif (A(m,n-1)==inf && A(m,n+1)==inf)
               e_x(m,n)=0;
           else
               e_x(m,n)=0;
           end
        end
        
        %Y-direction
        if (A(m,n)~=inf && m>1 && m<a)
           if (A(m-1,n))==inf
               e_y(m,n)=A(m,n)-A(m+1,n);
           elseif (A(m+1,n)==inf)
               e_y(m,n)=A(m-1,n)-A(m,n);
           else
               e_y(m,n)=(A(m-1,n)-A(m+1,n))/2;
           end
        end
            
        if (A(m,n)==inf && m>1 && m<a)
           if (A(m-1,n)==inf || A(m+1,n)==inf)
               e_y(m,n)=0;
           elseif (A(m-1,n)==inf && A(m+1,n)==inf)
               e_y(m,n)=0;
           else
               e_y(m,n)=0;
           end
        end
        
        if (n==1 || n==b)
            e_x(m,n)=0;
        end
        if (m==1 || m==a)
            e_y(m,n)=0;
        end
        
    end
end

        
       
        