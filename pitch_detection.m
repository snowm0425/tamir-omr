function text_pos = pitch_detection(diff_c,coord_crop,out_p,predicted,lThickness,lSpace)
%detect the pitch level of given symbol 'out_p'
%input:
%-diff_c: loaded cropdiff_%d.pbm
%-coord_crop: relative coordinate to the corresponding stave line
%-out_p: cropped symbol without staff lines
%-predicted: class prediction from svm
%-lThickness: thickness of staff line
%-lSpace: distance between two staff lines
%output:
%-text_pos: pitch level of input symbol
    
class_name = [cellstr('barline'),'clef_c','clef_f','clef_f_2','clef_g','color_breve','color_semibreve',...
'custos','fermata','flat','note_breve','note_fusa_down','note_fusa_up', ...
'note_longa_down','note_longa_up','note_maxima_down','note_maxima_up','note_minim_down','note_minim_up', ...
'note_semibreve','note_semiminim_down','note_semiminim_up','point','rest_breve','rest_longa','rest_minim', ...
'rest_semibreve','time_sig_Imin','time_sig_Imincut','time_sig_Min2','time_sig_Pmin','time_sig_Pmincut','time_sig_Triple'];
if(~isempty(strfind(class_name{predicted},'note')) && coord_crop(1,4)/coord_crop(1,3)>2) 
    %%note_fusa,longa,minim,semiminim
    text_mapping = cal_text_mapping(diff_c, coord_crop(1,1), lThickness,lSpace);
    [val,indx_line] =max(sum(out_p'));

     %detect bar position
    ysum_outp = sum(out_p');
    nonz_outp = find(ysum_outp>0);
    temp_outp = ysum_outp(nonz_outp(1):nonz_outp(end)); 
    temp_outp(temp_outp<mean(temp_outp)) = 0;
    if(indx_line>size(out_p,1)/2)  %up
        non_zero = find(temp_outp);
        start_head = non_zero(1);
        indx_nonz = 0;
        while(start_head<size(ysum_outp,2)/2)
            indx_nonz = indx_nonz+1;
            start_head = non_zero(indx_nonz);
        end
        middle_head = round((size(temp_outp,2)-start_head)/2)+start_head+nonz_outp(1);
    else
        non_zero = find(temp_outp);
        %     hal_non_zero = non_zero(1:round(size(non_zero,2)/2));
        start_head = non_zero(end); 
        indx_nonz = 0;
        while(start_head>size(ysum_outp,2)/2)
            indx_nonz = indx_nonz -1;
            start_head = non_zero(end+indx_nonz);
        end
        middle_head = round((start_head+nonz_outp(1))/2);  
    end
    if(round((indx_line+middle_head)/2)+coord_crop(1,2)<1)
        text_pos = text_mapping(1);
    else
        text_pos = text_mapping(round((indx_line+middle_head)/2)+coord_crop(1,2));
    end
elseif(~isempty(strfind(class_name{predicted},'note')) || ~isempty(strfind(class_name{predicted},'color')))
    %% note breve, colorred breve or semibreve
    text_mapping = cal_text_mapping(diff_c, coord_crop(1,1), lThickness,lSpace);
    val =(sum(out_p'));
    indx_line = round((find(val~=0,1)+find(val~=0,1,'last'))/2); %the middle point
    text_pos = text_mapping(indx_line+coord_crop(1,2)); 

elseif(predicted == 2) %clef c
    text_mapping = cal_text_mapping(diff_c, coord_crop(1,1), lThickness,lSpace);
    val =(sum(out_p'));
    firstOfMaxes = find(val>mean(val),1);
    lastsOfMaxes = find(val>mean(val),1,'last');
    %[maxValue, linearIndexesOfMaxes] = max(val(:));
    % [rowsOfMaxes firstOfMaxes] = find(val == maxValue,1,'first');
    % [rowsOfMaxes lastsOfMaxes] = find(val == maxValue,1,'last');
    indx_line = round((firstOfMaxes+lastsOfMaxes)/2);
    text_pos = text_mapping(indx_line+coord_crop(1,2)); 
elseif(predicted == 3) %clef f
    if(coord_crop(1,1)<=0)
        text_mapping = cal_text_mapping(diff_c,round(coord_crop(1,1)+coord_crop(1,3)/2), lThickness, lSpace);
    else
        text_mapping = cal_text_mapping(diff_c, coord_crop(1,1), lThickness,lSpace);
    end
    %select the right part
	[B,L] = bwboundaries(out_p,'noholes');
	mid = [];
	for idx_B=1:size(B,1)
        min_B = min(B{idx_B});
        if(min_B(2)>size(out_p,2)/2)
            max_B = max(B{idx_B});
            mid = [mid;round((max_B(1)+min_B(1))/2)];
        end
    end
	if isempty(mid) || numel(mid)==1
        text_pos = 100
    else
         %locate the middle of two dots
        middle_head = round((mid(1)+mid(2))/2); 
        text_pos = text_mapping(middle_head+coord_crop(1,2));
    end
elseif(predicted == 4) %clef f_2
	text_mapping = cal_text_mapping(diff_c, coord_crop(1,1), lThickness,lSpace);  
	%select the right part
	val =(sum(out_p));
    max_val = find(val==max(val),1);
	left_x = find(val(1:mv)==0,1,'last');
	right_val = sum(out_p(:,lef_x:size(out_p,2))');
  	%detect the bar line
	nonz_outp = find(right_val>0);
	temp_outp = right_val(nonz_outp(1):nonz_outp(end));
	temp_outp(temp_outp<mean(temp_outp)) = 0;
	non_zero = find(temp_outp);
	start_head = non_zero(end); 
	indx_nonz = 0;
	while(start_head>size(right_val,2)/2)
        indx_nonz = indx_nonz -1;
        start_head = non_zero(end+indx_nonz);
    end
	%locate the middle point
	middle_head = round((start_head+nonz_outp(1))/2); 
	text_pos = text_mapping(round((indx_line+middle_head)/2)+coord_crop(1,2));
elseif(predicted == 5) %clef g
	text_mapping = cal_text_mapping(diff_c, coord_crop(1,1), lThickness,lSpace);
	val = sum(out_p');
	indx_line = find(val~=0,1,'last');
	text_pos = text_mapping(indx_line+coord_crop(1,2))+0.5;		
elseif(predicted == 10) %flat
	text_mapping = cal_text_mapping(diff_c, coord_crop(1,1), lThickness,lSpace);
	val =(sum(out_p'));
	 %         lastsOfMaxes = find(val>mean(val),1,'last');
	%        firstOfMaxes = find(val == val(lastsOfMaxes),1);
    %       indx_line = round((firstOfMaxes+lastsOfMaxes)/2);
	indx_line = find(val~=0,1,'last');  
	text_pos = text_mapping(indx_line+coord_crop(1,2))+0.5; 
else % for the symbol which does not need the pitch level
	text_pos = -10; 
end

end
