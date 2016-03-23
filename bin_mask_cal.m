function mask= bin_mask_cal(diff_c,lThickness)
% build the binary mask to avoid the detection of lyrics/other ornaments
% output a binary mask in the same size as diff_c
 
en = 0;
mask = zeros(size(diff_c));
ref_start = [];
for i=1:size(diff_c,2)
   [staff_start,staff_end] = check_five_lines(diff_c(:,i),lThickness);
   if(size(staff_start,1)==5)
       mask(staff_start(1):staff_start(end),i) = 1;
       if ~en
           ref_start = i;
       end
   elseif sum(diff_c(:,i))>=lThickness
       mask(:,i) = 1;
   end
end
if isempty(ref_start)
    disp('no five staff lines found...');
end
tmp_sum = sum(mask);
% start from the ref_start to correct

full_list = find(tmp_sum==size(mask,1));
ref_list = find(full_list>ref_start);

for i=1:length(ref_list)
   indx = full_list(ref_list(i));
   if(sum(mask(:,indx-1))<size(mask,1) && sum(mask(:,indx-1))>0)
      first = find(mask(:,indx-1),1);
      last = find(mask(:,indx-1),1,'last');
      nonz = find(diff_c(:,indx));
      if(abs(nonz(1)-first)<lThickness)
         mask(1:nonz(1)-1,indx) = 0 ; 
      else
          mask(1:first-1,indx) = 0;
      end
      
      if(abs(nonz(end)-last)<lThickness)
          mask(nonz(end)+1:size(mask,1),indx) = 0;
      else
          mask(last+1:size(mask,1),indx) = 0;
      end
   elseif sum(mask(:,indx-1))==0
       mask(:,indx) = 0;       
   end
end


% go back to the non-five lines case

full_list = find(tmp_sum==size(mask,1));
ref_list = find(full_list<ref_start);
for i=length(ref_list):-1:1
   indx = full_list(ref_list(i));
   if(sum(mask(:,indx+1))<size(mask,1) && sum(mask(:,indx+1)>0))
      first = find(mask(:,indx+1),1);
      last = find(mask(:,indx+1),1,'last');
      nonz = find(diff_c(:,indx));
      if(abs(nonz(1)-first)<lThickness)
         mask(1:nonz(1)-1,indx) = 0 ; 
      else
          mask(1:first-1,indx) = 0;
      end
      
      if(abs(nonz(end)-last)<lThickness)
          mask(nonz(end)+1:size(mask,1),indx) = 0;
      else
          mask(last+1:size(mask,1),indx) = 0;
      end
   elseif sum(mask(:,indx+1))==0
       mask(:,indx) = 0;
   end
end
