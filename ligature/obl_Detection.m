
function obliqPosition = obl_Detection(lineOutput, w,lSpace)
% detect the obliques from the output of LSD detector
%% input: lineOutput : text output from ./lsd with staff line removed input
%% output: linePosition: return vertical line segments in n by 3( x, y1, y2) matrix, 
    obliqPosition = []; 
    fileID = fopen(lineOutput,'r');
    tline = fgetl(fileID);
    linePosition = [];
    potential_lines = [];
    avg_thickness = 0;
    num_segments = 0;
    while(ischar(tline))
        num_segments = num_segments+1;
        nums = strsplit(tline);
     %  x_y = str2num(char(nums(1:5)));
       x_y = [min(str2num(char(nums(1))),str2num(char(nums(3))));...
           min(str2num(char(nums(2))),str2num(char(nums(4)))); ...
          max(str2num(char(nums(1))),str2num(char(nums(3))));...
          max(str2num(char(nums(2))),str2num(char(nums(4)))); ...
          round(str2num(char(nums(5))))]; %x1,y1,x2,y2 
       bar_length = sqrt((x_y(1)-x_y(3))^2+(x_y(2)-x_y(4))^2);
       bar_angle = atand(abs((x_y(2)-x_y(4))/(x_y(1)-x_y(3))));
       if(bar_angle>5 && bar_angle<80) %group the lines between 5-80degree
          potential_lines = [potential_lines  [bar_length; x_y]];%length,x,y,x,y,thick
       end
       %linePosition = [linePosition; x_y' ]; %#ok<AGROW>
       avg_thickness = avg_thickness+x_y(5);
       tline = fgetl(fileID);
    end
    if ~isempty(potential_lines)
        
        avg_thickness = avg_thickness/num_segments;
        [temp_y,temp_i] = sort(potential_lines(2,:)); %sort x 
        candidate_lines = [];
        %filter out the line segments thinner than avg
        for i=1:size(potential_lines,2)
           if potential_lines(6,i)>avg_thickness
               candidate_lines = [candidate_lines ,potential_lines(:,i)];
           end
        end

        %merge the line segments in the same x range
        line_dist = zeros(1,round(w));
       % dif = candidate_lines(:,3)-candidate_lines(:,1);
       % max_indx = find(dif==max(dif));
       % line_dist(candidate_lines(max_indx,1):candidate_lines(max_indx,3)) = 1;
        if~isempty(candidate_lines)

            for i=1:size(candidate_lines,2)
               line_dist(candidate_lines(2,i):candidate_lines(4,i)) = 1; 
            end
            [B,L] = bwboundaries(line_dist,'noholes');

            if(size(B,1)==1) % contains one oblique
                min_x = max(min(B{1}));
                max_x = max(max(B{1}));
                obliqPosition = [min_x,max_x];
                disp('found one obliq');
            elseif(size(B,1)>1)
                disp('maybe contains more than one obliq');
                for i=1:size(B,1)
                   min_x = max(min(B{i}));
                   max_x = max(max(B{i}));
                   if(max_x-min_x)>0.5*lSpace
                      obliqPosition=[obliqPosition;min_x,max_x];
                %       disp('found one obliq');
                   end
                   
                   % obliqPosition = [obliqPosition;min_x,max_x];
                end
            end
        else
            disp('no candidate lines->no oblique');
        end
    else
       disp('no lines with angle between 10-80 degree detected->no oblique'); 
    end
    
    
