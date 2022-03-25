%boundary_value_editor
fn ='../data/v1v2v3';
roi = load(fn);

while(1)
figure(102); clf;
plot(roi.anchorpos(:,1),roi.anchorpos(:,2))
[x, y, button]=ginput(1);
if button==3
    break
end


dis = dot(roi.anchorpos - [x y], roi.anchorpos - [x y], 2);

[mindis,id]= min(dis);

if mindis < 0.1
    hold on;
    plot(roi.anchorpos(id,1),roi.anchorpos(id,2),'ro')    
    [x, y, button]=ginput(1);    
    if button==3
        break
    end
    roi.anchorpos(id,:) = [x y];
    
end


end
%%
anchor = roi.anchor;
anchorpos = roi.anchorpos;
id2delete = roi.id2delete;
save(fn, 'anchor', 'anchorpos','id2delete');