function parsing_and_recog(img_name,img_annotation,vlfeat)
%parse the segmentation file in annotation format and recognize the symbols with related information

%input:
%- img_name: the filename of the image
%- img_annotation: the output folder 
%- vlfeat: the path to where VLFEAT installed

%img_name = 'NLsHerAB_72A_003v';
%img_annotation = '/esat/jabbah/yhuang/test/ISMIR/';
%vlfeat = '/esat/jabbah/yhuang/vlfeat-0.9.19/';

currentFolder = pwd

% setup vl_feat
run(fullfile(vlfeat,'toolbox','vl_setup.m'));
addpath(genpath('./ligature'));

class_name = [cellstr('barline'),'clef_c','clef_f','clef_f_2','clef_g','color_breve','color_semibreve',...
'custos','fermata','flat','note_breve','note_fusa_down','note_fusa_up', ...
'note_longa_down','note_longa_up','note_maxima_down','note_maxima_up','note_minim_down','note_minim_up', ...
'note_semibreve','note_semiminim_down','note_semiminim_up','point','rest_breve','rest_longa','rest_minim', ...
'rest_semibreve','time_sig_Imin','time_sig_Imincut','time_sig_Min2','time_sig_Pmin','time_sig_Pmincut','time_sig_Triple'];

Params.PCA_COMP = 64;
Params.CLUSTERS = 128;
SVM_MODEL_FILE = fullfile(currentFolder,'data',sprintf('SVM_LIBLINEAR_MODEL_150sample_16gmm_FV_33class_noprob_%d_%d.mat',Params.PCA_COMP,Params.CLUSTERS));
ModelFile = fullfile(currentFolder,'data',sprintf('vl_gmm_model_expnote_pca_comp_16sample_v0_%d_clusters_%d.mat',Params.PCA_COMP,Params.CLUSTERS));

load(ModelFile); %model for computing fisher vector representation
load(SVM_MODEL_FILE); %model for classfying symbols
% load the y position for each stave
load(fullfile(img_annotation,img_name,'stave_position.mat'));

% load correct annotation file
fileID = fopen(fullfile(img_annotation,img_name,sprintf('%s_seg_col_correct.annotation',img_name)),'r');

outID = fopen(fullfile(img_annotation,img_name,sprintf('%s_recog.annotation',img_name)),'w');
im = imread(fullfile(img_annotation,img_name,sprintf('%s_crop.jpg',img_name)));
if ndims(im) == 3
    gray = rgb2gray(im);
else
    gray = im;
end

%load the Thickness parameter
path_DPfile = fullfile(img_annotation,img_name,'DP_log');
		
[status result] = system(['tail -n 1 ' path_DPfile]); %extract the last line of DP_log


if(strfind(result,'Thickness'))
    midstr = strsplit(result);
    thickstr = strrep(char(midstr(1)),',',' ');
	thickstr = strsplit(thickstr,':');
	spacestr = strsplit(char(midstr(2)),':'); 
	lThickness = str2double(thickstr(2));
	lSpace = str2double(spacestr(2));
end


en_proc = 0;
en_label = 0;
en_fv = 0;
tline = fgetl(fileID);
fprintf(outID,'%s\n',tline);
last_x = 500;
channel = 0;
which_line = 1; %defaut: start from first line

%---parse the annotation file from auto_crop_single.m---%
while ischar(tline)
   if(tline)
      if(tline(1)=='o') %id detected
         id = tline;  
         tline = fgetl(fileID)
         fprintf(outID,'%s\n',tline);
      elseif(tline(1)=='b') %bbox
         value = tline(7:end);
         split_v = strsplit(value,',');
         
         tline = fgetl(fileID);
         fprintf(outID,'%s\n',tline);
         en_proc = 1;    
      
      else
		%  fprintf(outID,'%s\n',tline);
          tline = fgetl(fileID);
          fprintf(outID,'%s\n',tline);
      end
      
      if(en_proc == 1)
          %determine which line
          if last_x-str2num(char(split_v(1)))>5*lSpace
             
             tmp_diff = stave_y(:,3) - round(str2num(char(split_v(2)))+str2num(char(split_v(4)))/2);
             tmp_diff = round(str2num(char(split_v(2))))-stave_y(:,2);
             tmp_diff(tmp_diff<0) = inf;
             which_line = find(tmp_diff==min(tmp_diff));
                if(size(which_line,1)~=1)
                    tmp_diff = abs(round(str2num(char(split_v(2))))-stave_y(:,2));
                    which_line = find(tmp_diff==min(tmp_diff));
                end
          end
          % ligature
          if(str2num(char(split_v(3)))>2*lSpace && str2num(char(split_v(1)))>last_x) %ligature
              disp('ligature');
              [label,pitch] = ligature_detection(img_name, img_annotation,str2num(char(split_v(1))),...
                  str2num(char(split_v(2))),str2num(char(split_v(3))),str2num(char(split_v(4))));
              if ~isempty(label)
                  fprintf('predict: %s\n',label); 
                  text_pos = num2str(pitch(1))
                  for i=2:size(pitch,2)
                    text_pos = strcat(text_pos,'_',num2str(pitch(i)))
                  end
                  fprintf(outID,'predicted: %s\nchannel: %d\ntext_position: %s\n\n',label,...
                  channel,text_pos);
                  en_fv = 0;
                  en_proc = 0;
                  
              else
                  en_fv = 1;
              end
          else
              en_fv = 1;
          end
          if en_fv
              tic;
              % classify cropped symbol in fisher vector representation
              FV_test = fv_calculate(gray,split_v,model);
              FV_test = double(FV_test ./ norm(FV_test));
              sparseh = sparse(FV_test); 
              [predicted,accuracy,dec_values] = predict(0,sparseh,svm_model);
              fprintf('predicted: %s \n',class_name{predicted});
              toc
              en_proc = 0;

              
             % calculate relative coordinate
            
             %extract the symbol from the whole staff lines removed img
              imn = imread(fullfile(img_annotation,img_name,sprintf('%s_out.pbm',img_name)));
              diff_c = imread(fullfile(img_annotation,img_name,sprintf('cropdiff_%d.pbm',which_line)));
              coord_crop = [round(str2num(char(split_v(1))))-stave_y(which_line,1),...
                  round(str2num(char(split_v(2))))-stave_y(which_line,2),...
                  round(str2num(char(split_v(3)))),round(str2num(char(split_v(4))))];
              out_p = imcrop(imn,[round(str2num(char(split_v(1)))),round(str2num(char(split_v(2)))),...
                  round(str2num(char(split_v(3)))),round(str2num(char(split_v(4))))]);
            
              if(isempty(out_p)) disp('empty out_p'); end
              
              if(coord_crop(1,1)<=0)
                  coord_crop(1,1) = round((coord_crop(1,1)+coord_crop(1,4))/2);
              end
             % if coord_crop(1,2)<0
             %     coord_crop(1,2)=1
             % end

             %%-------- pitch detection---------%
              
             text_pos = pitch_detection(diff_c,coord_crop,out_p,predicted,lThickness,lSpace);
              fprintf(outID,'predicted: %s\nchannel: %d\ntext_position: %f\n\n',class_name{predicted},channel,text_pos);
              if predicted == 1
                    channel = channel+1; 
              end
          end
          last_x = str2num(char(split_v(1)));
          
      end
   else
       tline = fgetl(fileID);
       fprintf(outID,'%s\n',tline);
   end
end
fclose(outID);
fclose(fileID);
end
       
