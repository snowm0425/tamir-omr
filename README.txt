TAMIR - OMR

00-Introduction

TAMIR is a toolbox for transcriping music manuscript from 16th century into modern music format. It mainly contains two parts, OMR and transcription. This package contains the OMR part of the toolbox. The transcription part can be found in 
https://github.com/zb8c5ek/alamire_of_tamir


01-Dependencies

-Matlab

-VLFEAT(http://www.vlfeat.org/install-matlab.html)

-LIBLINEAR(https://www.csie.ntu.edu.tw/~cjlin/liblinear/)

-Piotr's Computer Vision Matlab Toolbox(http://vision.ucsd.edu/dollar/toolbox/doc/index.html)

-LSD line segment detector(http://www.ipol.im/pub/art/2012/gjmr-lsd/)
	Note: compile the source code and put the executable(lsd) under ligature folder

For visualization:
-Image Annotation Tool(https://github.com/alexklaeser/imgAnnotation)

For compiling our source code (in /StaffRemoval):
-qmake

02-Package structure

The package tamir/ is organized as follows:
- README.txt	: this file
- LICENSE 	: GPL v3 license
- data/ 	: a directory containing pre-trained models
- img/ 		: a directory containing an example input image (NLsHerAB_72A_003v.jpg)
- StaffRemoval/ : a directory containing source code of staff lines removal
- training/ 	: a directory containing Matlab scripts for training
- ligature/	: a directory containing source codes for ligature recognition
- output/	: a directory containing an example corrected segmentation file under /NLsHerAB_72A_003v folder for demo

03-Installation 
- move the excutable(lsd) form LSD line segment detector under /tamir
- move the compiled mex file of LIBLINEAR for your system architecture (ex:64bit, 'predict.mexa64') under /tamir
- go into /StaffRemoval and do the following:
	$ qmake -o Makefile DP_proc.pro
	$ make
It will generate an executable called 'DP_StaffRemoval'

04-License

This toolbox is distributed under the GPL v3 license (see LICENSE.txt). If you use
this code for research purposes, please cite our paper [1] .

05-Reference
[1]TAMIR Toolbox: Recognition and Transcription Tools for Handwritten Mensural Manuscripts (under reviewing)

06-Quick start
Install the libraris in 01-Dependencies and follow 03-Installation steps.

------Manuscript preprocessing + symbol segmentation-------
To run the example for preprocessing and symbol segmentation, start Matlab and type
>> preprocessing_segmentation 'NLsHerAB_72A_003v' /path/to/tamir/img /path/to/tamir/output
Output file('NLsHerAB_72A_003v_seg_col.annotation') can be found in 
'/path/to/tamir/output/NLsHerAB_72A_003v/' directory.

------Symbol recognition--------
For symbol recognition, type
>> parsing_and_recog 'NLsHerAB_72A_003v' /path/to/tamir/output '/path/to/<VLFEATROOT>/' '/path/to/<Piotr>/toolbox/' '1'
It will load pre-corrected segmentation file and output the recognition result ('NLsHerAB_72A_003v_correct_recog.annotation') in
'/path/to/tamir/output/NLsHerAB_72A_003v/' directory.

------Visualization of the segmentation/recognition result-----------
Go into the directory where Image Annotation Tool is installed, start the software by
$ ./imgAnnotation

Specify the annotation file by "Database" -> "Open Database"
and select the annotation file to display.

Click the + sign in File/Directory column and click the image file, the segmentation/recognition result together with the image source will show up.

07-Contact
Will be disclosed after paper reviewing.

