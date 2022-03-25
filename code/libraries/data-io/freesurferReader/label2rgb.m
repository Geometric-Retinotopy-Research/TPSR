function color=label2rgb(label)

bin=dec2bin(label,24);

color(:,1)=bin2dec(bin(:,1:8));

color(:,2)=bin2dec(bin(:,9:16));


color(:,3)=bin2dec(bin(:,17:24));

color=double(color)/255;

end