function classifier_training
%script for training the symbol classifier

%%%%%%%%%%%%%USER DEFINE%%%%%%%%%%%%%%%%%%%%%%%%%%
%VLFEAT set up; remember to change the path to vl_setup.m
run('/path/to/vlfeat-0.9.19/toolbox/vl_setup.m');
%specify the location of the training images
alldata_dir = '/path/to/training/';
%specify the location of images to train a GMM
Params.image_dir = '/path/to/GMM_imgs16/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Params.data_dir = './data/';

%------------------------------%
fv_phow = 1;

labels = [];
training = 1;
if(training)
	if(fv_phow)
		% parameter
		Params.PCA_COMP = 64;
		Params.CLUSTERS = 128;
		Params.phowOptsNeg = {'Step', 3} ;
		Params.keyPointsPerImage = 300;
		Params.sift_samples4gmm = 10000;
		SVM_MODEL_FILE = fullfile(Params.data_dir,sprintf('SVM_LIBLINEAR_MODEL_150sample_16gmm_FV_33class_noprob_%d_%d.mat',Params.PCA_COMP,Params.CLUSTERS));

		% ------------------------ create gmm model -------------------------------
		model = createGMMModel(Params);
		%-------------------------- create fisher vectors -------------------------
		fprintf('***** Start processing fisher vectors ****\n');
		if exist(SVM_MODEL_FILE,'file') ~= 2
			FISH_FILE = fullfile(Params.data_dir,sprintf('Fisher_Vec_33class_150sample_16gmm_v0_%d_%d_4SVMLIN.mat',Params.PCA_COMP,Params.CLUSTERS));
			if exist(FISH_FILE,'file') == 2
				load(FISH_FILE);
      
			else
				class_dir = dir(alldata_dir);
				FV_Pos = [];
				label = [];
				for i= 3: size(class_dir,1)
					data_dir = fullfile(alldata_dir,class_dir(i).name,'unprocessed');
					temp_FV_Pos = createFisher(data_dir,model,Params.phowOptsNeg);
					FV_Pos = [FV_Pos ; temp_FV_Pos];
					size_training = size(temp_FV_Pos,1);
					label = [label ; zeros(size_training,1)+i-2];
				end
				save(FISH_FILE,'FV_Pos','label','-v7.3');%}
			end
        
			FV_Pos_sparse = sparse(FV_Pos);

    	    %--------------------------- train SVM or what ever -----------------------
			fprintf('***** Start training  ****\n');
			svm_model = train(label, FV_Pos_sparse);
			save(SVM_MODEL_FILE,'svm_model','-v7.3');
		else
			load(SVM_MODEL_FILE);
		end
	end
end

end


%--------------------functions--------------------------------------%

function model = createGMMModel(Params)
	ModelFile = fullfile(Params.data_dir,sprintf('vl_gmm_model_expnote_pca_comp_16sample_v0_%d_clusters_%d.mat',Params.PCA_COMP,Params.CLUSTERS));
	SIFT_Data_File = fullfile(Params.data_dir,sprintf('sift_data_expnote_16sample_file_v0_%d_%d.mat',Params.keyPointsPerImage,Params.sift_samples4gmm));
	if exist(ModelFile,'file') ~= 2
		if exist(SIFT_Data_File,'file') ~= 2
			SIFT = getSiftForGMMCreation(Params);
			save(SIFT_Data_File,'SIFT','-v7.3');
		else
			load(SIFT_Data_File);
		end
		[pcamap] = princomp(single(SIFT'));
		SIFT =  pcamap(:,1:Params.PCA_COMP)' * single(SIFT);
		[means, sigmas, weights] = vl_gmm(SIFT, Params.CLUSTERS);
		model.pcamap = pcamap;
		model.means = means;
		model.sigmas = sigmas;
		model.weights = weights;
		save(ModelFile,'model','-v7.3');
     
	else
		load(ModelFile);
	end
 
end

function SIFT_data = getSiftForGMMCreation(Params)
	ims = dir(fullfile(Params.image_dir, '*.jpg'))' ;
	num_files = size(ims,2);       
	fprintf('num_files:%d',num_files);
	index = 1;    
	kpPerIm = Params.keyPointsPerImage;
	sift_samples = Params.sift_samples4gmm;
	MAX_Sam = num_files * kpPerIm;
	fprintf('sift_samples:%d',sift_samples);
	SIFT_data = zeros(128,MAX_Sam,'uint8');
	for i = 1 : num_files
		im = imread(fullfile(Params.image_dir,ims(i).name)) ; 
		if ndims(im) == 3
			im = rgb2gray(im);
		end  
		[drop, siftrn] = vl_phow(single(im), Params.phowOptsNeg{:}) ;        
		siftrn = siftrn(:,sum(siftrn) ~= 0 ) ;
		if size(siftrn,2) > kpPerIm
			rn = randperm(size(siftrn,2));
			siftrn = siftrn(:,rn(1:kpPerIm));
		end
		SIFT_data(:,index:(index+size(siftrn,2)-1)) = siftrn;
		index = index + size(siftrn,2);
		fprintf('index:%d\n',index);
		if(index ==  MAX_Sam+1)
			fprintf('Exceeded %d,index:%d,MAX=%d\n',i,index,MAX_Sam);
			break
		end
		fprintf('Progress %d  / %d \n',i,num_files);
	end    
	rn = randperm(size(SIFT_data,2));
	SIFT_data = SIFT_data(:,rn(1:sift_samples));
end


%--------------------------------------------------------------------------
% Creates fisher vectors for a given input image directory dirName from randomly selected 150 %images determined by permutation_d.bin in each directory
% model is the GMM model and pcamap
% phowOpts is the vl_phow options
% return column fisher vectors
function FV = createFisher(dirName,model,phowOpts)
%--------------------------------------------------------------------------
	PCA_COMP = size(model.means,1);
	CLUSTERS = size(model.means,2);
	ims = dir(fullfile(dirName, '*.jpg'))' ;
	count_t = 0;
	num_files = size(ims,2);
	if(num_files>150)
		num_files = 150;
	end
	FV = zeros(num_files, PCA_COMP * CLUSTERS * 2 ,'double');
	fprintf('size of FV:');
	disp(size(FV));
	SiftSampleSize = 10000;
	tt = tic();
	bin_p = fullfile(dirName,'permutation_d.bin'); %random permutation file for each class
	fileID = fopen(bin_p,'r');
	permu_id = fread(fileID,'double');
	fclose(fileID);
	for i = 1 : size(ims,2)
		t = tic();
		count_t = count_t+1;
		
		if permu_id(i)<=num_files
			file = fullfile(dirName,ims(i).name);
			count_t = count_t+1;
		else
			continue;
		end
		
		im = imread(file) ; 
		im = im2single(im) ;
		if ndims(im) == 3
			im = rgb2gray(im);
		end 
		if(size(im,1)>size(im,2))
			length = size(im,1); 
		else
			length = size(im,2);
		end
		im = imresize(im, 240/length,'bilinear'); 
		im_h = size(im,1);
		im_w = size(im,2);
		[drop, siftrn] = vl_phow(single(im), phowOpts{:}) ; 
		remove_ix = sum(siftrn) ~= 0 ;
		drop = drop(:, remove_ix) ;
		siftrn = siftrn(:, remove_ix ) ;
		siftrn = model.pcamap(:,1:PCA_COMP)' * single(siftrn);
		ENC = vl_fisher(siftrn	, model.means, model.sigmas, model.weights, 'improved');
		temp_ENC = ENC';
		fprintf('size of tmpENC');
		disp(size(temp_ENC));
		FV(count_t,:) = temp_ENC;
		FV(count_t,:) = FV(count_t,:) ./ norm(FV(count_t,:));
		FV(count_t,:) = double(FV(count_t,:));
		t = toc(t);
		ttt = toc(tt);
		fprintf('%d :%d elapsed : %1.2f   Total : %1.2f......%d estimated : %1.2f \n',count_t,i,t,ttt,num_files,(ttt/i)*(num_files-i));
	end
	if(count_t == num_files)
		disp('equal!');
	else
		disp('not equal;');
	end
end

