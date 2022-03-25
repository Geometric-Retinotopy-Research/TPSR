function y=predictbold( pp, stimdd, hrf, res, xx,yy)
% Specifically, the variable <pp> is a vector of parameter values (1 x 5)
% and the variable <dd> is a matrix with the stimuli (frames x pixels).

resmx = max(res);
% pre-compute some cache

global keyframe keyframe2full ddsim  

if(isempty(keyframe))
    fprintf('\n Pre-Calculate the mapping from keyframe 2 full frames');
    try 
        load('keyframe_keyframe2full.mat','keyframe', 'keyframe2full','ddsim' );
    catch
   
    for i = 1:size(stimdd,1) 
        foundsame = 0 ;
        for j = length(keyframe):-1:1
            
            framei = stimdd(i,:);
            framej = stimdd(keyframe(j),:);
            if ( quick_isequal(  framei, framej)) % the purpose is quick check for unequal                
               
                
                foundsame =1;
                break  % if equal, which means it is a duplicated one
                
            end
            
        end 
        if(foundsame)
            keyframe2full(i) =  j;
        else
            keyframe = [keyframe; i]; 
            keyframe2full(i) =  length(keyframe);
        end
        
        if (mod(20*i, size(stimdd,1) ) ==0)
            fprintf('.');
        end
    end
        
    ddsim = stimdd(keyframe,:);
    
    save('keyframe_keyframe2full.mat','keyframe', 'keyframe2full','ddsim' );
    end
end
 
 
% stdgauss = makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0); %Elapsed time is 0.002887 seconds.
stdgauss = makegaussian2d_sef(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0); %Elapsed time is 0.002887 seconds.

placedstd = stdgauss/(2*pi*abs(pp(3))^2);

expsth = ddsim*[vflatten(placedstd); 0]; %Elapsed time is 0.006948 seconds.

rawy = posrect(pp(4)) * expsth.^ posrect(pp(5)); % Elapsed time is 0.000767 seconds.

rawyfull = rawy(keyframe2full); % Elapsed time is 0.000030 seconds.
 
y = conv2run_sef(rawyfull, hrf); % Elapsed time is 0.002283 seconds.

end


function f=conv2run_sef(y, hrf)

f = zeros([1800,1],'single');
for i=0:300:1800-1
    temp = conv2(y(i+1:i+300),hrf);
    f(i+1:i+300) = temp(1:300,:);
end
 
end



function flag = quick_isequal( fi, fj)

flag = sum(fi~=fj)<10;

end
 


function f = makegaussian2d_sef(res,r,c,sr,sc,xx,yy,ang,omitexp)  
r = (-1/res) * r + (.5 + .5/res);  % this is faster
c = (1/res) * c + (-.5 - .5/res);  % this is faster
sc = sc/res;
 f = exp(((xx-c).^2+(yy-r).^2)/(-(2*sc^2)));
end