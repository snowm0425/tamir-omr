function text_mapping = cal_text_mapping(diff_c,coord_tempx,lThickness,lSpace)
%obtain the staff lines mapping from a given x position for pitch level detection 
add_value = 0;
uconv = [1,1,1,1,1];
[staff_start,staff_end] = check_five_lines(diff_c(:,(coord_tempx)),lThickness);
while size(staff_start,1)~= 5			 	
    if(add_value >=0)
		if(add_value <size(diff_c,2)-coord_tempx-5*lSpace)
		   add_value = add_value + 1;
        else
           add_value = -1;
		end
	else                
		add_value = add_value - 1;
	end
	[staff_start,staff_end] = check_five_lines(diff_c(:,coord_tempx+add_value),lThickness);
end
%replace convolution
line_sum = zeros(size(diff_c,1),1);
for iter_staff = 1:5
    diff_lines = staff_end(iter_staff)-staff_start(iter_staff);
    line_sum((staff_start(iter_staff)-diff_lines):(staff_end(iter_staff)+diff_lines)) = 1; 
end
divided_thickness = round(max(line_sum)/5);
%diff_conv = conv(uconv,single(diff_c(:,coord_tempx+add_value)));
%diff_conv(diff_conv~=0)=1;
%diff_conv = diff_conv(3:end-2);
list_one = find(line_sum==1);
if(list_one(end)+round(lSpace/2)>size(line_sum,1))
    above_one = single(line_sum(1:list_one(end)+round(size(line_sum,1)-list_one(end))/2));
else
    above_one = single(line_sum(1:list_one(end)+round(lSpace/2)));
end
above_one(above_one==0) = 0.5;
above_one(above_one==1) = 0;
text_mapping = [above_one;zeros((size(diff_c,1)-size(above_one,1)),1)];
for ind_abv = (size(above_one,1)-1):-1:1
   if above_one(ind_abv+1)~=above_one(ind_abv)
	  text_mapping(ind_abv) = text_mapping(ind_abv+1) + 0.5;
   else
      text_mapping(ind_abv) = text_mapping(ind_abv+1);
   end
end 		

