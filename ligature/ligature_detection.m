
function [label,pitch] = ligature_detection(img_name, img_annotation,coord_x, coord_y, coord_w, coord_h, ptoolbox)
% recognize the ligature from a given cropped symbol
% xx im: crop of ligature(rgb)
% coord_x, coord_y, coord_w, coord_h : bbox for the symbol
% ptoolbox: path to where Piotr's toolbox installed (/path/to/piotr/toolbox/)
% output label: ex: Ligature_uB_B_B
%        pitch: ex: 3.5_4.5_3


img = imread(fullfile(img_annotation,img_name,sprintf('%s_crop.jpg',img_name)));
im = imcrop(img,[coord_x,coord_y,coord_w,coord_h]);
path_DPfile = fullfile(img_annotation,img_name,'DP_log');
entire_w = zeros(1,round(coord_w));		
[status result] = system(['tail -n 1 ' path_DPfile]); %extract the last line of DP_log


if(strfind(result,'Thickness'))
    midstr = strsplit(result);
    thickstr = strrep(char(midstr(1)),',',' ');
	thickstr = strsplit(thickstr,':');
	spacestr = strsplit(char(midstr(2)),':'); 
	lThickness = str2double(thickstr(2));
	lSpace = str2double(spacestr(2));
end

% block detection
addpath(genpath(ptoolbox))

% load pre-trained modelfor note breve det
load('./data/breved5Detector.mat');
%detector = acfTrain( opts );
tic, bbs=acfDetect(im,detector); toc
new_bbs = [];
sorted_bbs = [];
if(~isempty(bbs))
    num_breve = round(size(im,2)/mean(bbs(:,3)));
    %make all component positive
    negative_ind = find(bbs<0);
    if ~isempty(negative_ind)
       for i = 1:length(negative_ind) 
        bbs(negative_ind(i)) = 1;
       end
    end
    if(size(bbs,1)>num_breve) 
       %start candidate selection
       tmp_indx = find(bbs(:,5)==max(bbs(:,5)));
       new_bbs = bbs(tmp_indx,:);
       bbs(tmp_indx,:)=[]; %delete the row with max score
       num_out = 1;
       %%%%%TO DO!!!!!  %%%%
       while(~isempty(bbs) && num_out<num_breve)
           tmp_indx = find(bbs(:,5)==max(bbs(:,5)));
           tmp_cor = bbs(tmp_indx,:);
           if(tmp_cor(1,5)>90) %score >90
              %determin whether to take, check the overlapping ratio with current candidates
              en_take = check_overlap(new_bbs,tmp_cor);
              if(en_take)
                 new_bbs = [new_bbs;tmp_cor]; 
              end
           end
           bbs(tmp_indx,:)=[];
       end
        
    end
    if(size(new_bbs,1)>1)
    	[d1,d2] = sort(new_bbs(:,1))
        sorted_bbs = new_bbs(d2,:);
        % bounding box refinement!!
         %TODO

         %append the symbol name, pitch
        sorted_bbs = [sorted_bbs,ones(size(sorted_bbs,1),1)*'B',ones(size(sorted_bbs,1),1)*100]; %later check for color

        for i=1:size(sorted_bbs,1)
            entire_w(round(sorted_bbs(i,1)):round(sorted_bbs(i,1)+sorted_bbs(i,3))) = 1; 
        end
    else
        sorted_bbs = [];
    end
end


[B,W]=bwboundaries(~entire_w,'noholes');
check_obl = 0;
if ~isempty(B)
    for i =1:size(B,1)
        min_x = max(min(B{i}));
        max_x = max(max(B{i}));
        if(max_x-min_x > lSpace)
            check_obl = 1;
        end
    end
end
if check_obl
    img_bw = imread(fullfile(img_annotation,img_name,sprintf('%s_out.pbm',img_name)));
    im_bw = imcrop(img_bw,[coord_x,coord_y,coord_w,coord_h]);
    imwrite(im_bw,'tmp.pgm');
    system('./lsd -s 0.5 tmp.pgm tbw.txt.result');
    system('rm -rf tmp.pgm');
    obliqXPosition = obl_Detection('tbw.txt.result',coord_w,lSpace);
    %check the obl detection result and return y positions for starting and
    %ending point
    % one oblique:[x1,y1,0,h1;x2,y2,0,h2]
    
    obliq_bbx = obl_y_Detection(obliqXPosition,im_bw,lSpace,lThickness);
    
    if isempty(obliq_bbx)
        check_obl = 0;
    elseif(size(obliq_bbx,1)>2)
        %use the height to eliminate the repetative obliq
        obliq_bbx = obl_repcheck(obliq_bbx,lThickness);
        
    end
    
    
end

% detect maxima?
% TODO


%pitch level calculation
pitch = [];

 % determine which line of stave it comes from
load(fullfile(img_annotation,img_name,'stave_position.mat'));
tmp_diff = stave_y(:,3) - coord_y;
tmp_diff(tmp_diff<0) = inf;
which_line = find(tmp_diff==min(tmp_diff));

diff_c = imread(fullfile(img_annotation,img_name,sprintf('cropdiff_%d.pbm',which_line)));

 if ~isempty(sorted_bbs)
    for num_breve = 1:size(sorted_bbs,1)
       text_mapping = cal_text_mapping(diff_c,...
           round(sorted_bbs(num_breve,1)+sorted_bbs(num_breve,3)/2+coord_x-stave_y(which_line,1)), lThickness,lSpace);
       indx_line = round(sorted_bbs(num_breve,2)+round(sorted_bbs(num_breve,4)/2));
       text_pos = text_mapping(round(indx_line+coord_y-stave_y(which_line,2)));
       % check the pitch
       indx_upper = text_mapping(round(sorted_bbs(num_breve,2)+coord_y-stave_y(which_line,2)));
       indx_lower = text_mapping(round(sorted_bbs(num_breve,2)+sorted_bbs(num_breve,4)+coord_y-stave_y(which_line,2)));
       if indx_upper-indx_lower==1
           pitch = [pitch;(indx_upper+indx_lower)/2];
           sorted_bbs(num_breve,7) = (indx_upper+indx_lower)/2;
       else
          pitch = [pitch;text_pos];
          sorted_bbs(num_breve,7) = text_pos;
       end
    end
 end

%pitch calculation for obl

if check_obl
  %combine the obl position to obliq_bbx with score, label: 'O', and pitch
  obliq_bbx = [obliq_bbx, zeros(size(obliq_bbx,1),1),ones(size(obliq_bbx,1),1)*'O', ones(size(obliq_bbx,1),1)*100];
  
  for num_obl = 1:size(obliq_bbx,1)
     %%% TODO if(obliq_bbx(num_obl,2) -->check obl_y_Detection!
      text_mapping = cal_text_mapping(diff_c,...
       round(obliq_bbx(num_obl,1)+coord_x-stave_y(which_line,1)), lThickness,lSpace);
   indx_line = round(obliq_bbx(num_obl,2)+round(obliq_bbx(num_obl,4)/2));
   text_pos = text_mapping(round(indx_line+coord_y-stave_y(which_line,2)));
   obliq_bbx(num_obl,7) = text_pos
  end
  % combine the result of breve and oblique and sort
  all_bbs = [sorted_bbs;obliq_bbx];
  [d1,d2] = sort(all_bbs(:,1))
    sorted_bbs = all_bbs(d2,:);
    %remove the wrong detected breve in between two 'O's
    obli_indx = find(sorted_bbs(:,6)=='O');
    if(obli_indx)
        for obli_i =1:2:size(obli_indx,1)
           if(obli_indx(obli_i)<size(sorted_bbs,1))
               %check the next label
               if(sorted_bbs(obli_indx(obli_i)+1,6)=='B')
                   %remove this line
                   sorted_bbs(obli_indx(obli_i)+1,:) = [];
               end
           end
        end
    end
end


% vertical line detection

% line detection %
if isempty(sorted_bbs)
    label=[];
    pitch = [];
else
    imwrite(im,'tmp.pgm');
    system('./lsd -s 0.5 tmp.pgm tmp.txt.result');
    system('rm -rf tmp.pgm');

    linePosition =  bar_Detection('tmp.txt.result',size(im,1));

    % determine the vertical lines or u or d or nothing
      %collect the starting x to a row matrix
      if(sorted_bbs(:,6) == 'B')
          out_position = [sorted_bbs(:,1)',sorted_bbs(end,1)+sorted_bbs(end,3)];
      else
          out_position = sorted_bbs(:,1)';
      end
      out_p = ones(1,size(sorted_bbs,1)*2+1)*'_';
      for i=1:size(sorted_bbs,1)
          out_p(2*i) = sorted_bbs(i,6);
      end
    %out_p = repmat('_',size(out_position));
    % create a matrix with '_' in between each character label for output label
    % replace '_' with 'u' or 'd'
    % eliminate '_' between two OO and in 1st and last position
    for i=1:size(linePosition,1)
       tmp = abs(out_position-linePosition(i,1)); 
       [d1 j] = min(tmp);
       if(j==1) %left most
           if(sorted_bbs(j,2)-linePosition(i,2)>size(im,1)/3)
              out_p(j) = 'u'
           elseif(linePosition(i,3)-sorted_bbs(j,2)-sorted_bbs(j,4)>size(im,1)/3)
               out_p(j) = 'd'
           end
       elseif(j==size(out_position,2)) %right most
           if(sorted_bbs(j-1,2)-linePosition(i,2)>size(im,1)/3)
               out_p(j) = 'u'
           elseif(linePosition(i,3)-sorted_bbs(j-1,2)-sorted_bbs(j-1,4)>size(im,1)/3)
               out_p(j) = 'd'
           end    
       else %middle
           if(min(sorted_bbs(j,2),sorted_bbs(j-1,2))-linePosition(i,2)>size(im,1)/3)
               out_p(j) = 'u'
           elseif(linePosition(i,3)-max(sorted_bbs(j,2)+sorted_bbs(j,4),sorted_bbs(j-1,2)+sorted_bbs(j-1,4))>size(im,1)/3)
               out_p(j) = 'd'
           end
       end

    end
    %label = char(out_p);
    %eliminate '_' between two 'O's
    obl_indx = find(out_p=='O');
    if obl_indx
        for i = 1:2:size(obl_indx,1)
            if obl_indx(i)<size(out_p,2)
               out_p(obl_indx(i)+1)= ' '; 
            end
        end
    end
    label = strrep(out_p, ' ','');
    if(label(end)=='_')
       label(end) = []; 
    end
    if(label(1)=='_')
        label(1) = [];
    end
    % return output label
    pitch = sorted_bbs(:,7)';
    label = strcat('Ligature_',label)
    disp(pitch);
end
    
    
    

end


function en= check_overlap(new_bbs, bbs)
    en = 0;
    for i = 1:size(new_bbs,1)
        if(new_bbs(i,1)>bbs(1,1))
            if(new_bbs(i,1)>bbs(1,1)+bbs(1,3))
                en = en+1;
            elseif((bbs(1,1)+bbs(1,3)-new_bbs(i,1))/new_bbs(i,3)<0.4 )
                en = en+1;
            end
        else
            if(new_bbs(i,1)+new_bbs(i,3)<bbs(1,1))
                en = en+1;
            else
                if((new_bbs(i,1)+new_bbs(i,3)-bbs(1,1))/new_bbs(i,3)<0.4 && bbs(1,1)+bbs(1,3)>new_bbs(i,1)+new_bbs(i,3))
                    en = en+1;
                end
            end
        end
        
    end
    if(en~=size(new_bbs,1))
        en=0; end

end
