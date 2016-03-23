function symbol_separate(path_DPfile, pbmin, pbmout,img_annotation,img_name)
%segment the symbols from the image using connected component analysis

[status result] = system(['tail -n 1 ' path_DPfile]); %extract the last line of DP_log


if(strfind(result,'Thickness'))
    midstr = strsplit(result);
	thickstr = strrep(char(midstr(1)),',',' ');
	thickstr = strsplit(thickstr,':');
	spacestr = strsplit(char(midstr(2)),':'); 
	lThickness = str2double(thickstr(2));
	lSpace = str2double(spacestr(2));

	%-------symbol_separate.....-----------%

	im_in = imread(pbmin);
	im_out = imread(pbmout);
	imdiff = im_in - im_out;
	se = strel('line',lThickness,0);
	%imdiff = imdilate(imdiff,se);
	flip = imdiff(:,round(size(imdiff,2)/2):size(imdiff,2));
	sum_y = sum(flip'); 
	y_axis = size(sum_y,2);
	sum_y(sum_y<=mean(sum_y)+50) = 0;
	sum_y(sum_y~=0) = 1;

	indx = find(sum_y,1);
	if sum_y(indx+lSpace+1) ==0 %%%
	   indx = indx+lSpace+find(sum_y(indx+lSpace:size(sum_y,2)),1);
	end
	staff = [indx];
	flag = 1;
	while(indx<y_axis)
	   indx = indx + lThickness*5+lSpace*6;
	   indx = indx + find(sum_y(indx:y_axis),1);
	   staff = [staff indx];
	end
	pattern = (lThickness+lSpace)*8
	stave = [];
	for i=1:size(staff,2)
	   disp(i);
	   crop = imcrop(imdiff,[0,staff(i)-2*(lThickness+lSpace),size(imdiff,2),pattern]);
	   sum_x = sum(crop);
	   sum_x(sum_x<lThickness*4)=0;
	   diff = find(sum_x,20);
	   cons_x = 0 
	   for j = 1:size(diff,2)-1
			disp(j); 
			if(diff(j+1)-diff(j)<3)
				cons_x = cons_x +1;
				disp('+1');
			else
				cons_x = 0;
			end
			if cons_x == 1
				indx_x = diff(j);
			elseif cons_x>2
				stave = [stave indx_x];
				break;
			end
	 
	   end
	end
	out = [];
	mask_stave = [];
		%img = imread('A-Wn_15495_002v_out.pbm');

	if(size(stave,2)>0 && size(stave,2)==size(staff,2))
	   for i=1:size(stave,2)
		  out = [out; (stave(i)-15), staff(i)-2*(lThickness+lSpace), size(imdiff,2),pattern];
		  mask_stave = [mask_stave;stave(i)-(lThickness+lSpace),staff(i)-lThickness/2,staff(i)+lThickness*3+4*(lThickness+lSpace)]; %x1,y1,y2
	  %refine width of the stave
		  diff_c = imcrop(imdiff,[(stave(i)-15),staff(i)-2*(lThickness+lSpace),size(imdiff,2),pattern]);
		  diffc_sum = sum(diff_c);
		  diffc_sum(diffc_sum<2*lThickness) = 0;
		  find_array = find(diffc_sum);
		  c = imcrop(im_out,[(stave(i)-15),staff(i)-2*(lThickness+lSpace),find_array(end),pattern]);
		%% end of width refinement  
		%% calculating mask for the stave %%
		  crop_h = size(diff_c,1);
		  crop_w = size(diff_c,2);
		  sum_crop_y = sum(diff_c);
		  non_zero_crop = find(sum_crop_y~=0);
		  en_finding_start = 1;
		  cons_crop = 0;
		  for q = non_zero_crop(1):non_zero_crop(end)
			  if(sum_crop_y(q)~=0)
                 if(en_finding_start == 1)
					start_crop_i = q;
					en_finding_start = 0;
                 else
    				cons_crop =  cons_crop+1;
	 			 end
              else
                 if(cons_crop>crop_w/2)
					end_crop_i = q;
                	break;
                 else	
                    cons_crop = 0;
                    en_finding_start = 1;
                 end
              end
		 end	
	%	 fprintf('start:%d,end:%d',start_crop_i,end_crop_i);
		 lin_indx = find(sum_crop_y==max(sum_crop_y),1); %%%%%%
		 lin_y = diff_c(:,lin_indx);
		 start_y_indx = find(lin_y~=0,1);
		 end_y_indx = find(lin_y~=0,1,'last');
		% binmask = zeros(size(diff_c));
        binmask = bin_mask_cal(diff_c,lThickness);
		%  if(start_y_indx>lThickness*2)
		 %     start_y_indx = start_y_indx-lThickness*2
		 % end
		 %binmask(start_y_indx:end_y_indx,start_crop_i:find_array(end))=1;
		 imwrite(binmask,fullfile(img_annotation,img_name,sprintf('mask_%d.pbm',i)));
		%% end of mask calculation
		 imwrite(c,fullfile(img_annotation,img_name,sprintf('cropin_%d.pbm',i)));
		 imwrite(diff_c,fullfile(img_annotation,img_name,sprintf('cropdiff_%d.pbm',i)))
       end
    end
	disp(out);
    % save the starting/ending y coordiate of each stave to a mat
    stave_y = [out(:,1),out(:,2),out(:,2)+out(:,4)]; % x , y , y2
    save(fullfile(img_annotation,img_name,'stave_position.mat'),'stave_y','-v7.3');
    
	if(isempty(out))
        %fprintf(fileLog,'%s:staff extraction failed.\n',img_name);
		disp('staff extraction failed!')
	else

	%----- extract symbols ----%
		fileAno = fopen(fullfile(img_annotation,img_name,sprintf('%s_seg_col.annotation',img_name)),'w');
		fprintf(fileAno,'########## NEW FILE ##########\nfile: %s_crop.jpg\n\n',fullfile(img_annotation,img_name,img_name));
		fileAnox = fopen(fullfile(img_annotation,img_name,sprintf('%s_seg2_col.annotation',img_name)),'w');
		fprintf(fileAnox,'########## NEW FILE ##########\nfile: %s_crop.jpg\n\n',fullfile(img_annotation,img_name,img_name));
		coord_diffx = out(:,1,:);
		coord_diffy = out(:,2,:);
		global_i = 0;
		global_j = 0;
		for num_crop = 1:size(coord_diffx,1)
			last_predicted = 0;
			
			imn = imread(fullfile(img_annotation,img_name,sprintf('cropin_%d.pbm',num_crop)));
			diff_c = imread(fullfile(img_annotation,img_name,sprintf('cropdiff_%d.pbm',num_crop)));
			mask = imread(fullfile(img_annotation,img_name,sprintf('mask_%d.pbm',num_crop)));
			se = strel('line',lThickness,45);
			%im = imdilate(imn,se);
			im =imn;
			[B,L] = bwboundaries(im,'noholes');
			coord = [];
			coordy = [];
			coord_xy = [];
			coord_xy_imo = [];
			mid_coord = [];
			text_pos = 0;
			add_value = 0;
			last_x1 = 0;
			last_y1 = Inf; %to avoid the decoration included
			last_x2 = 0;
			last_y2 = 0;
			up_en = 0;
			down_en = 0;
			last_top = 0;
			last_bottom = 0;
			for i = 1:length(B)
				max_B = max(B{i});
				min_B = min(B{i});
			   	if(max_B(2)-min_B(2)>lThickness && max_B(1)-min_B(1)>lThickness && length(B)>1 && min_B(2)+coord_diffx(num_crop)<size(im_in,2)-lSpace)
					if(sum(sum(mask(min_B(1):max_B(1),min_B(2):max_B(2))))>0) %inside the mask %%%%%%%
						global_j = global_j+1;	
						coord_xy_imo = [coord_xy_imo;coord_diffx(num_crop)+min_B(2),coord_diffy(num_crop)+min_B(1),max_B(2)-min_B(2),max_B(1)-min_B(1)];	
						fprintf(fileAnox,'object: %d\nbbox: %d,%d,%d,%d\nposition:in\n\n',global_j,coord_diffx(num_crop)+min_B(2),coord_diffy(num_crop)+min_B(1),max_B(2)-min_B(2),max_B(1)-min_B(1));
						last_y1 = min_B(1);
						last_y2 = max_B(1);
                        last_x1 = min_B(2);
                        last_x2 = max_B(2);
						
					elseif(max_B(1)<size(imn,1)/2 && max_B(1)>last_y1) % above the stave
						global_j = global_j+1;	
						coord_xy_imo = [coord_xy_imo;coord_diffx(num_crop)+min_B(2),coord_diffy(num_crop)+min_B(1),max_B(2)-min_B(2),max_B(1)-min_B(1)];	
						fprintf(fileAnox,'object: %d\nbbox: %d,%d,%d,%d\nposition:above\n\n',global_j,coord_diffx(num_crop)+min_B(2),coord_diffy(num_crop)+min_B(1),max_B(2)-min_B(2),max_B(1)-min_B(1));
						last_y1 = min_B(1);
						last_y2 = max_B(1);
                        last_x1 = min_B(2);
                        last_x2 = max_B(2);
					elseif(max_B(1)>size(imn,2)/2 && min_B(1)<last_y2 && max_B(1)<last_y2)
						global_j = global_j+1;	
						coord_xy_imo = [coord_xy_imo;coord_diffx(num_crop)+min_B(2),coord_diffy(num_crop)+min_B(1),max_B(2)-min_B(2),max_B(1)-min_B(1)];	
						fprintf(fileAnox,'object: %d\nbbox: %d,%d,%d,%d\nposition:below\n\n',global_j,coord_diffx(num_crop)+min_B(2),coord_diffy(num_crop)+min_B(1),max_B(2)-min_B(2),max_B(1)-min_B(1));
						last_y1 = min_B(1);
						last_y2 = max_B(1);
                        last_x1 = min_B(2);
                        last_x2 = max_B(2);
					end			
				else	
					continue;
			   	end	
				
            end
            
			if ~isempty(coord_xy_imo)
				coord_temp = [coord_xy_imo(1,:)];
				for i = 2:size(coord_xy_imo,1)
					last_y = coord_xy_imo(i-1,2);
					last_x = coord_xy_imo(i-1,1);
					last_h = coord_xy_imo(i-1,4);
					last_w = coord_xy_imo(i-1,3);
					last_mid_x = last_x+last_w/2;
					last_mid_y = last_y+last_h/2;
					cur_mid_x = coord_xy_imo(i,1)+coord_xy_imo(i,3)/2;
					cur_mid_y = coord_xy_imo(i,2)+coord_xy_imo(i,4);
					if(last_w/coord_xy_imo(i,3)<1.45 && last_w/coord_xy_imo(i,3)>0.6 && coord_xy_imo(i,1)<=last_x+last_w+lThickness/2 && cur_mid_y>last_y && cur_mid_y-last_mid_y < lSpace)	
						if(max(last_x+last_w,coord_xy_imo(i,1)+coord_xy_imo(i,3))-min(last_x,coord_xy_imo(i,1))>2*lSpace) % avoid two closed notes 
							coord_temp = [coord_temp;coord_xy_imo(i,:)];
						else
							coord_temp(end,:) = [min(last_x,coord_xy_imo(i,1)),min(coord_xy_imo(i,2),last_y),max(last_x+last_w,coord_xy_imo(i,1)+coord_xy_imo(i,3))-min(last_x,coord_xy_imo(i,1)),max(last_y+last_h,coord_xy_imo(i,2)+coord_xy_imo(i,4))-min(coord_xy_imo(i,2),last_y)];
                            coord_xy_imo(i,:) = coord_temp(end,:);
						end
					elseif(i==size(coord_xy_imo,1) && coord_xy_imo(i,4)>4*lSpace && last_h>4*lSpace && coord_xy_imo(i,3)<lSpace/2 && last_w<lSpace/2) %barline
						coord_temp(end,:) = [last_x,last_y,coord_xy_imo(i,3)+coord_xy_imo(i,1)-last_x,coord_xy_imo(i,2)+coord_xy_imo(i,4)-last_y];
					elseif(i==size(coord_xy_imo,1) && coord_xy_imo(i,1)-last_x-last_w>3*lSpace)
						fprintf('noise on edge!\n');
					else
						coord_temp = [coord_temp;coord_xy_imo(i,:)];
					end
                end	
				
                %merge the oversegmented boxes along the same x-axis
                coord_temp2 = coord_temp(1,:);
                for indxx=2:length(coord_temp)
                    if(coord_temp(indxx,1)<coord_temp2(end,1)+coord_temp2(end,3)) &&...
                        (coord_temp2(end,1)+coord_temp2(end,3)-coord_temp(indxx,1))/coord_temp(indxx,3)>0.6 &&...
                        coord_temp(indxx,2)<coord_temp2(end,2)+coord_temp2(end,4)     
                            coord_temp2(end,:) = [min(coord_temp(indxx,1),coord_temp2(end,1)),...
                                min(coord_temp(indxx,2),coord_temp2(end,2)),...
                                max(coord_temp(indxx,1)+coord_temp(indxx,3),coord_temp2(end,1)+coord_temp2(end,3))-...
                                min(coord_temp(indxx,1),coord_temp2(end,1)),...
                                max(coord_temp(indxx,2)+coord_temp(indxx,4),coord_temp2(end,2)+coord_temp2(end,4))-...
                                min(coord_temp(indxx,2),coord_temp2(end,2))];
                    else
                            coord_temp2 = [coord_temp2;coord_temp(indxx,:)];
                    end
                end
                
                
                
				for i = 1:size(coord_temp2,1)
					global_i = global_i+1;
					coord_temp2(i,1) = coord_temp2(i,1)-lThickness;
					coord_temp2(i,2) = coord_temp2(i,2)-lThickness;
					coord_temp2(i,3) = coord_temp2(i,3)+2*lThickness;
					coord_temp2(i,4) = coord_temp2(i,4)+2*lThickness;
					%%%% pitch detection %%%
					fprintf(fileAno,'object: %d\nbbox: %d,%d,%d,%d\n\n',global_i,coord_temp2(i,:));
				end		
			end
			
		end
		fclose(fileAno);
		fclose(fileAnox);
    end
else
    %fprintf(fileLog,'DP proc. failed.\n');
	disp('staff removal fails!')
end
