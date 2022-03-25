%% read a subject path !!!!!!!!!!!!!!!!!!!! Note: Not every files is readed
%% Syntax
% subject = read_subject_59k_fs_LR(subjectpath)

%% Description
% Inputs:
%       1. subjectpath: subject path
% Outputs:
%       1. subject: struct

%% Contribution
%  Author : Yanshuai Tu
%  Created: 2019/xx/xx
%  Revised: Not yet
%
%  Copyright @ Geometry Systems Laboratory
%  School of Computing, Informatics, and Decision Systems Engineering, ASU
%  http://gsl.lab.asu.edu/

%% Reference :

function subject = read_subject_59k_fs_LR(subjectpath)
% 
fprintf('Reading subject from path %s .', subjectpath);

% analysis subjectid
files = dir(sprintf('%s/*.gii',subjectpath));
giisamplefn = files(1).name;
dotids=strfind(giisamplefn,'.');
subjectid = giisamplefn(1:dotids(1)-1);

% read  surfaces
feature_names= {'inflated', 'sphere', 'pial', 'white', 'midthickness', 'very_inflated', ...
    'thickness', 'curvature', 'aparc', 'sulc', 'BA'};
feature_types ={'surf','shape','label'};

for i =1:length(feature_names)
    
    errcnt =0;
    for j= 1: length(feature_types)
        try
            g=gifti(sprintf('%s/%s.L.%s.59k_fs_LR.%s.gii', subjectpath, subjectid, feature_names{i}, feature_types{j}));
            if (isfield(g,'cdata'))
                g= g.cdata;
            end
            eval(sprintf('subject.lh.%s = g;', feature_names{i}));
            
            g=gifti(sprintf('%s/%s.R.%s.59k_fs_LR.%s.gii', subjectpath, subjectid, feature_names{i}, feature_types{j}));
            if (isfield(g,'cdata'))
                g= g.cdata;
            end
            eval(sprintf('subject.rh.%s = g;', feature_names{i}));
            break;
        catch
            errcnt = errcnt +1;
            
            if(errcnt ==length(feature_types))
                error(sprintf('Read error for feature %s ',feature_names{i}))
            end
        end
    end
    fprintf('.');
end
 
fprintf('Finish Reading\n');

end