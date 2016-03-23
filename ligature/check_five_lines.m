function [staff_start,staff_end] = check_five_lines(diff_c_line, lThickness)
%check the input column vector contains five lines or not
uconv = [1,1,1,1,1];
%{
[Gmag,Gdir] = imgradient(diff_c_line);
mean_g = mean(Gmag(Gmag>0));
std_g = std(Gmag(Gmag>0));
Gmag(Gmag<(mean_g-std_g)) = 0;
Gmag(Gmag~=0) = 1;
tmp_line = conv(uconv,Gmag);
%build a mask to avoid noise
mask_three = tmp_line;
mask_three(mask_three~=3)=0;
list_three = find(mask_three);
if ~isempty(list_three)
   for ind_t = 1:length(list_three)
       if list_three(ind_t)-lThickness>0
          mask_three(list_three(ind_t)-lThickness:list_three(ind_t)+lThickness) = 1;
       elseif list_three(ind_t)+lThickness > length(tmp_line)
           mask_three(list_three(ind_t)-lThickness: length(tmp_line)) = 1;
       else
           mask_three(1:list_three(ind_t)) = 1;
       end
   end
end

tmp_line(tmp_line~=0) = 1;
tmp_line = tmp_line .* mask_three;
line_sum = tmp_line;
%}

staff_start = [];
staff_end = [];
cum_thick = 1;
list_one = find(diff_c_line);
if(~isempty(list_one))
    tmp_start = list_one(1);
    tmp_end = list_one(1);
    for ind_o = 2:length(list_one)
        if(list_one(ind_o)-tmp_start<lThickness*1.5 && ind_o~=length(list_one))
            cum_thick = cum_thick+1;
            tmp_end = list_one(ind_o);
        else
            if(abs(cum_thick-lThickness)<lThickness/2)
                cum_thick = 1;
                staff_end = [staff_end;tmp_end];
                staff_start = [staff_start;tmp_start];
            end
            tmp_start = list_one(ind_o);
        end
    
    end
    
end

