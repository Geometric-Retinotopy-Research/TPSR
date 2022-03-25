clear all;
fprintf('Compiling mex files ... ');

% compile mex files
mex mex/perform_front_propagation_anisotropic.cpp

disp('done.');