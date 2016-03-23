function FV_test = fv_calulate(gray,split_v,model)
%calculate the fisher vector representation for given intpu image 'gray'
%input: 
%-gray: gray scaled full image,
%-split_v: coordinate of cropped symbol
%-model: pre-trained GMM model
%output:
%-FV_test: fisher vector representation for given cropped symbol
    
Params.phowOptsNeg = {'Step', 3} ;

im_crop = gray(round(str2num(char(split_v(2)))):round(str2num(char(split_v(2)))+str2num(char(split_v(4)))),...
round(str2num(char(split_v(1)))):round(str2num(char(split_v(1)))+str2num(char(split_v(3)))));
en_label = 0;
if(size(im_crop,1)>size(im_crop,2))
	length = size(im_crop,1); 
else
	length = size(im_crop,2);
end
			  
im_crop = imresize(im_crop, 240/length,'bilinear'); %fix longest length to 240
im_h = size(im_crop,1);
PCA_COMP = size(model.means,1);
CLUSTERS = size(model.means,2);
phowOpts = Params.phowOptsNeg;

FV_test = zeros(1, PCA_COMP * CLUSTERS * 2 ,'double');
[drop, siftrn] = vl_phow(single(im_crop),phowOpts{:}) ; 
remove_ix = sum(siftrn) ~= 0 ;
drop = drop(:, remove_ix) ;
siftrn = siftrn(:, remove_ix ) ;
siftrn =  model.pcamap(:,1:PCA_COMP)' * single(siftrn);

ENC = vl_fisher(siftrn, model.means, model.sigmas, model.weights, 'improved');
FV_test = ENC';
clear ENC;

end
