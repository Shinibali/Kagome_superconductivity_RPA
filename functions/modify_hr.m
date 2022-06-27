function [] = modify_hr()

fwan = './input_files/'; 
fname = 'hoppings_w90style.txt'; fout = 'modified_hr.mat';

fid = fopen([fwan,fname]);
lines = textscan(fid,'%s','delimiter','\n'); fclose(fid); 
lines = lines{1}; nol = max(size(lines));

Rvec = zeros(nol-2,2) ; sisj = zeros(nol-2,2) ;
oioj = zeros(nol-2,2) ; HR = zeros(nol-2,1);
for i = 3:nol
    temp = regexp(strtrim(lines{i}),'\s+','split'); %strsplit(strtrim(lines{i}),' ');
    if cellfun(@isempty,temp) ~= 1 
        Rvec(i-2,1:2) = [str2double(temp{1}), str2double(temp{2})] ; 
        sisj(i-2,1:2) = [str2double(temp{3}), str2double(temp{4})] ;
        oioj(i-2,1:2) = [str2double(temp{5}), str2double(temp{6})] ;
        HR(i-2) = str2double(temp{7}) ; 
    end
end
norb = max(max(oioj)); nsite = max(max(sisj));

inputs = [nol-2, 7]; frc = 'NrowNcol_tbm.bin';
fid = fopen([fwan,frc],'w+','n');fwrite(fid,inputs,'double'); fclose(fid);
inputs = [Rvec,sisj,oioj,HR]; fhop = 'hopping.bin';
fid = fopen([fwan,fhop],'w+','n');fwrite(fid,inputs,'double'); fclose(fid);

save([fwan,fout],'norb','nsite','Rvec','sisj','oioj','HR','-mat')
fclose('all');

end