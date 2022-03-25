function [fx_bound fy_bound]=boundary(A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for computation of boundary forces.         %
% A-input matrix (obtained from provagetfile.m)        %
% fx_bound - boundary force in x-direction             %
% fy_bound - boundary force in y-direction             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parameters
   %Parameter of boundary potential
   u0_alphaB=10;
   R=0.2;

%A=A(11:85,3:133);
A=flipud(A);
[a b]=size(A);
wall=find(A==5);
A1=ones(a,b);
A1(wall)=0;
[y x]=find(A1==0);

for m=1:a
    for n=1:b
        if (A1(m,n)~=0)
            rx=n-x;
            ry=m-y;
            rr=sqrt(rx.^2+ry.^2);
            r(m,n)=min(rr);
        end
        if (A1(m,n)==0)
            r(m,n)=0;   
        end
    end
end

r(r==0)=-inf;
bound_pot=u0_alphaB.*exp(-r./R);
[fx_bound fy_bound]=grad(bound_pot);
figure,quiver(fx_bound,fy_bound,2)

figure
quiver(fx_bound,fy_bound,2),hold on
[y_wall x_wall]=ind2sub(size(A1),wall)
for i=1:size(x_wall)
    plot(x_wall(i),y_wall(i),'.k','LineWidth',1.5),hold on
end



