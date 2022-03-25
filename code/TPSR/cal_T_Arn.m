function [T,Arn] = cal_T_Arn(Froi, u, v, V)

bd = compute_bd(Froi);
fa = face_area(Froi, u);
fv = face_area(Froi,v);
fu = face_area(Froi,u);
Froi(ismember(Froi(:,1),bd)| ismember(Froi(:,2),bd) |ismember(Froi(:,3),bd) | fv<1e-4| fu<1e-4,: )=[];
flip_fid = find(abs(compute_bc(Froi, u,v))>1);
T = length(flip_fid);
Arn = sum(face_area(Froi(flip_fid,:), V)) / sum(face_area(Froi, V) );
end