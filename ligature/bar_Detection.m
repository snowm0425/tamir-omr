
function linePosition = bar_Detection(lineOutput, h)
%detect the vertical lines from the LSD detector output

%% input: lineOutput : text output from ./lsd with staff line removed input
%% output: linePosition: return vertical line segments in n by 3( x, y1, y2) matrix, 
%%         where n is the number of vertical lines 

    fileID = fopen(lineOutput,'r');
    tline = fgetl(fileID);
    linePosition = [];
    vertical_lines = [];
    while(ischar(tline))
       nums = strsplit(tline);
     %  x_y = str2num(char(nums(1:5)));
       x_y = [min(str2num(char(nums(1))),str2num(char(nums(3))));...
           min(str2num(char(nums(2))),str2num(char(nums(4)))); ...
          max(str2num(char(nums(1))),str2num(char(nums(3))));...
          max(str2num(char(nums(2))),str2num(char(nums(4)))); ...
          round(str2num(char(nums(5))))]; %x1,y1,x2,y2 
       bar_length = sqrt((x_y(1)-x_y(3))^2+(x_y(2)-x_y(4))^2);
       bar_angle = atand(abs((x_y(2)-x_y(4))/(x_y(1)-x_y(3))));
       if(abs(bar_angle-90)<10) %group the lines from 90 degree
          vertical_lines = [vertical_lines  [bar_length; x_y]];
       end
       %linePosition = [linePosition; x_y' ]; %#ok<AGROW>
       tline = fgetl(fileID);
    end
    [temp_y,temp_i] = sort(vertical_lines(2,:));
    avg_thickness = mean(vertical_lines(6,:))
    sorted_lines = vertical_lines(:,temp_i);
    if(size(sorted_lines,2)>0)
       tmp_array = zeros(round(max(max(sorted_lines(3,:)),max(sorted_lines(5,:)))),...
           round(max(max(sorted_lines(2,:)),max(sorted_lines(4,:)))));
        for i = 1:size(sorted_lines,2)
            if((sorted_lines(2,i)-sorted_lines(6,i))>0 && (sorted_lines(4,i)+sorted_lines(6,i))<size(tmp_array,2))
                 tmp_array(sorted_lines(3,i):sorted_lines(5,i),...
                    (sorted_lines(2,i)-round(sorted_lines(6,i)/2)):(sorted_lines(4,i)+round(sorted_lines(6,i)/2)))=1; 
            end
        end
    end
     [B,L] = bwboundaries(tmp_array,'noholes');
     temp_B = [];
     for i = 1:length(B)
        max_B = max(B{i});
        min_B = min(B{i});
     %   if(max_B(1)-min_B(1)>h/3)
            temp_B=[temp_B;round((min_B(2)+max_B(2))/2),min_B(1),max_B(1)]; % extend for img with stave
      %  end
     end 
     out_index = 1;
     linePosition = temp_B(1,:);
     for i = 2:size(temp_B,1)
         if(temp_B(i,1)>linePosition(out_index,1)+avg_thickness)
             linePosition = [linePosition;temp_B(i,:)];
             out_index = out_index+1;
         else
             if(temp_B(i,2)<linePosition(out_index,2))
                 linePosition(out_index,2) = temp_B(i,2);
             end
             if(temp_B(i,3)>linePosition(out_index,3))
                 linePosition(out_index,3) = temp_B(i,3);
             end
         end
         
     end
     
    %linePosition = tmp_array;
end
