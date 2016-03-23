function preprocessing_segmentation(img_name,img_source,img_annotation)
% preprocessing, binarization, staff lines removel, symbol segmentation

%-----------parameters------------%
%img_name = 'NLsHerAB_72A_003v';
currentFolder = pwd
%img_source = fullfile(currentFolder,'img');
%img_annotation = fullfile(currentFolder,'data'); %output directory
resize_ratio = 1;  % percentage for resize, 1 means the original size
%color_en = 0;
crop = 1;
#fileLog = fopen('./log','w');

%---------------------------------%
if(~exist(img_annotation,'dir'))
	mkdir(img_annotation);
end
if(~exist(fullfile(img_annotation,img_name),'dir'))
	mkdir(fullfile(img_annotation,img_name));
end

if(~exist(fullfile(img_annotation,img_name,sprintf('%s_seg.annotation',img_name)),'file'))
 % no segmentation file
    imgr_name = img_name
	%---cropping the margin and resize---%
    if crop
        if ~exist(fullfile(img_annotation,img_name,sprintf('%s_crop.jpg',img_name)),'file')
            crop_im = im_crop_xy(img_source,img_name,img_annotation,resize_ratio);
        else
            crop_im = imread(fullfile(img_annotation,img_name,sprintf('%s_crop.jpg',img_name)));
        end
    else
        crop_im = imread(fullfile(img_annotation,img_name,sprintf('%s.jpg',img_name)));
		if resize_ratio<1
			crop_im = imresize(crop_im,resize_ratio);
        end
    end
	%-----binarization-------%

    disp('binarizing...');
	level = graythresh(crop_im);
	bw = im2bw(crop_im,level);
	o = ~bw;
	disp('binarized');

	%-----color-------%
%    if color_en ==1
%       % if(~exist(fullfile(img_annotation,img_name,sprintf('%s.pbm',imgr_name)),'file'))
%       if(~exist(fullfile(img_annotation,img_name,'color.pbm'),'file'))
%            disp('color mask calculating..');
%	        fill = color_mask_cal(crop_im);
%            fill = ~fill;
%			imwrite(fill,fullfile(img_annotation,img_name,sprintf('color.pbm')),'pbm');
%       else
%            fill = imread(fullfile(img_annotation,img_name,'color.pbm'));
%       end   
%       if(sum(sum(fill))/size(fill,1)/size(fill,2)>0.7) % in case of color filtering fails
%            o = o .* fill;
%       end	
%    end
    imwrite(o,fullfile(img_annotation,img_name,sprintf('%s.pbm',imgr_name)),'pbm');
	
	%------staff lines removal-----%
	pbmin = fullfile(img_annotation,img_name,sprintf('%s.pbm',imgr_name));
	pbmout = fullfile(img_annotation,img_name,sprintf('%s_out.pbm',imgr_name));
	path_DPfile = fullfile(img_annotation,img_name,'DP_log');
	if(~exist(pbmout,'file'))
		disp('Staff lines removing ...');
	
		system(['./StaffRemoval/DP_proc ' pbmin ' ' pbmout ' >' path_DPfile]);
    end

	%-----symbol segmentation-----%
    symbol_separate(path_DPfile,pbmin,pbmout,img_annotation,img_name);
	
		
end
%fclose(fileLog);

end
