function crop_im = im_crop_xy(img_source,img_name,img_annotation,resize_ratio)

%remove the non-score part, ex: the ruler or the background when scanning the manuscript
%output: cropped image


%----resize the image---%
if resize_ratio<1
   im = imread(fullfile(img_source,sprintf('%s.jpg',img_name)));
   im = imresize(im,resize_ratio);
   imgr_name = sprintf('%s_%dp',img_name,resize_ratio*100);
   imwrite(im,fullfile(img_annotation,img_name,sprintf('%s.jpg',imgr_name)));
else
   im = imread(fullfile(img_source,sprintf('%s.jpg',img_name)));
  
   imwrite(im,fullfile(img_annotation,img_name,sprintf('%s.jpg',img_name)));
end
        
gray = rgb2gray(im);
bw = edge(gray);
sum_bw = sum(bw);
sum_bw(sum_bw<10) = 0;
accum = 0;
start_indx = 0;
end_indx = 0;
for i = 2:size(sum_bw,2)
	if(accum > 300)
		if sum_bw(i)~=0
			accum = accum+1;
		else
  			end_indx = i;
    		continue;
		end  
	else
		if(sum_bw(i-1)>0) && sum_bw(i)>0
            if accum == 0
                start_indx = i-1;
			end
			accum = accum +1;
        else
            accum = 0;
		end
    end
end

fprintf('start: %d, end: %d',start_indx,end_indx);

crop_im = imcrop(im,[start_indx-10,1,end_indx-start_indx,size(im,1)]);
%------crop horizontally(x-axis)-------%
gray = rgb2gray(crop_im);
m = mean(gray);
list_m = find(m<100,1);
if(isempty(list_m))
	%fprintf(fileLog,'%s:crop failed.\n',img_name);
	display('crop failed!')
elseif(list_m<size(gray,2)/2) %crop left part
	crop_im = imcrop(crop_im,[list_m,1,size(gray,2)-list_m,size(gray,1)]);
else
	crop_im = imcrop(crop_im,[1,1,list_m,size(gray,1)]);
end
%------crop in y-axis--------%
gray = rgb2gray(crop_im);
m=mean(gray');
start_y = find(m>100,1);
end_y = find(m>100,1,'last');
crop_im = imcrop(crop_im,[1,start_y,size(gray,2),end_y-start_y]);
imwrite(crop_im,fullfile(img_annotation,img_name,sprintf('%s_crop.jpg',img_name)));
disp('image cropped!');
