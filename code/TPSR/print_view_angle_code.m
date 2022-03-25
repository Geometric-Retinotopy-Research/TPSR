%%
% Description  -- function print_view_angle_code()
%       print the information of the current view
% 
%%
function print_view_angle_code()
va=camva;
pos=get(gca,'CameraPosition');

fprintf('camva(%f);\n',va)

fprintf('set(gca,''CameraPosition'',[%f %f %f])\n',pos(1),pos(2), pos(3))

pos=get(gca,'Position');
if length(pos) ==3
    fprintf('set(gca,''Position'',[%f %f %f %f])\n',pos(1),pos(2), pos(3))
else
    fprintf('set(gca,''Position'',[%f %f %f %f])\n',pos(1),pos(2), pos(3), pos(4))
end

b = get(gca,'Xlim');
c = get(gca,'Ylim');
d = get(gca,'Zlim');
fprintf('set(gca,''Xlim'',[%f %f])\n',b(1), b(2))
fprintf('set(gca,''Ylim'',[%f %f])\n',c(1), c(2))
fprintf('set(gca,''Zlim'',[%f %f])\n',d(1), d(2))

pos=get(gcf,'position');


fprintf('set(gcf,''position'',[%f %f %f %f])\n',pos(1),pos(2), pos(3), pos(4))


hLegend = findobj(gcf, 'Type', 'Legend');
pos=get(hLegend,'position');
try
fprintf('hLegend = findobj(gcf, ''Type'', ''Legend'');\n');

fprintf('set(hLegend,''position'',[%f %f %f %f])\n',pos(1),pos(2), pos(3),pos(4))
catch
    
end

end