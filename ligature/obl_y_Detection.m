function obl_bbx = obl_y_Detection(obliqXPosition,im_bw,lSpace,lThickness)
%check the y-axis of the candidate obliques from obliqXPosition 
%input: 
%output: obl_bbx: ex: one oblique:[x1,y1,0,h1;x2,y2,0,h2]
obl_bbx = [];
%%TOCHECK!! Inf component
for i = 1:size(obliqXPosition,1) % search for all x candiates
    %check the x starting point of obliq
    x1_slice = im_bw(:,obliqXPosition(i,1));
    y1=find(x1_slice,1);
    y2=find(x1_slice,1,'last');
    
    while(y2-y1>1.5*lSpace)
        x1_slice = im_bw(:,obliqXPosition(i,1)+lThickness);
        y1=find(x1_slice,1);
        y2=find(x1_slice,1,'last');
        obliqXPosition(i,1) = obliqXPosition(i,1)+lThickness;
        if obliqXPosition(i,1)+lThickness>=size(im_bw,2)
           y1=[];
           y2=[];
        end
    end
    j=0;
    p1=y1;
    p2=y2;
  %  p1=find(im_bw(:,obliqXPosition(i,1)-j),1);
 %   p2=find(im_bw(:,obliqXPosition(i,1)-j),1,'last');
   if(isempty(y1) || isempty(y2))
       obl_bbx = [];
   else
     while(abs(p1-y1)<lThickness || abs(p2-y2)<lThickness)
            j=j+1;
            y1=p1;
            y2=p2;
            if(obliqXPosition(i,1)>j)
                p1=find(im_bw(:,obliqXPosition(i,1)-j),1);
                p2=find(im_bw(:,obliqXPosition(i,1)-j),1,'last');
                if isempty(p1) || isempty(p2) || y2-y1>2*lSpace || abs(p1-y1)>lThickness || abs(p2-y2)>lThickness
                    p1=Inf;
                    p2=Inf;
                    y1 = find(im_bw(:,obliqXPosition(i,1)-j+2),1);
                    y2 = find(im_bw(:,obliqXPosition(i,1)-j+2),1,'last');
                    if isempty(y1) || isempty(y2)
                       y1 = 0;
                       y2 = 0;
                    end
                end
            else
                j=j-1;
                p1=Inf;
                p2=Inf;
            end
        end
        x1 = obliqXPosition(i,1)-j+lThickness;
         %check the x ending point
        x2_slice = im_bw(:,obliqXPosition(i,2));
        y3=find(x2_slice,1);
        y4=find(x2_slice,1,'last');
        while(y2-y1>1.5*lSpace)
            x1_slice = im_bw(:,obliqXPosition(i,1)-lThickness);
            y1=find(x1_slice,1);
            y2=find(x1_slice,1,'last');
            obliqXPosition(i,1) = obliqXPosition(i,1)-lThickness;
        end
        j=0;
        p1=y3;
        p2=y4;
        %p1=find(im_bw(:,obliqXPosition(i,2)+j),1);
        %p2=find(im_bw(:,obliqXPosition(i,2)+j),1,'last');
        if(isempty(y3) || isempty(y4))
            obl_bbx=[];
        else
            while(abs(p1-y3)<lThickness || abs(p2-y4)<lThickness)
                j=j+1;
                y3=p1;
                y4=p2;
                p1=find(im_bw(:,obliqXPosition(i,2)+j),1);
                p2=find(im_bw(:,obliqXPosition(i,2)+j),1,'last');
                 if isempty(p1) || isempty(p2) || y4-y3>2*lSpace || abs(p1-y3)>lThickness || abs(p2-y4)>lThickness
                    p2=Inf;
                    p1=Inf;
                    y3= find(im_bw(:,obliqXPosition(i,2)+j-2),1);
                    y4 = find(im_bw(:,obliqXPosition(i,2)+j-2),1,'last');
                end
            end
            x2 = obliqXPosition(i,2)+j-lThickness;
            obl_bbx = [obl_bbx;x1,y1,0,y2-y1;x2,y3,0,y4-y3];
        end
   end
end
%check/delete the obliq with small height
if ~isempty(obl_bbx)
    new_bbx = [];
   for i = 1:2:size(obl_bbx,1) 
       if(obl_bbx(i,4)>0.75*lSpace && obl_bbx(i+1,4)>0.75*lSpace)
           new_bbx = [new_bbx;obl_bbx(i,:);obl_bbx(i+1,:)];
       end
   end
  obl_bbx = new_bbx;
   
end
