function [fx_dest,fy_dest]=destination(e_x,e_y,v_x,v_y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used for determination of destination   %
% forces. Force has a form of an acceleration term and for %
% its determination one needs components of desired dire-  %
% ction vector (e_x,e_y) which are obtained using FMA and  %
% actual velocity vector (v_x,v_y).                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parameters
  %relaxation time
     tau_alpha=0.5;
  %mean desired velocity
     v0_mean=1.34;

%normalization of desired direction
norm=sqrt(e_x.^2+e_y.^2);
e_x=e_x./norm;
e_y=e_y./norm

deviation1=find(e_x==inf);
deviation2=find(e_y==NaN);

e_x(deviation1)=0;
e_x(deviation2)=0;

e_y(deviation1)=0;
e_y(deviation2)=0;

%destination forces
fx_dest=(1/tau_alpha).*(v0_mean.*e_x-v_x);
fy_dest=(1/tau_alpha).*(v0_mean.*e_y-v_y);

