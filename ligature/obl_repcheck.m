
function new_bbx =  obl_repcheck(obliq_bbx,lThickness)
%eliminate repetative candidates of oblique

for i = 1:2:size(obliq_bbx,1)
    current_h = obliq_bbx(i,4);
    if(abs(obliq_bbx(i+1,4)-current_h)>5*lThickness)
        obliq_bbx(i,:) = 0;
        obliq_bbx(i+1,:) = 0;
    end
    
end
ind=find(obliq_bbx(:,1)==0);
obliq_bbx(ind,:) = [];
new_bbx = obliq_bbx;
