%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to deal with FS contour post-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = FScontour_postproc(klist,bandind,len_dkf,fdir,fout,b1b2)

rearrklist = klist; rearrbdin = bandind; len_element = len_dkf ;
len_element = [len_element, zeros(max(size(len_element)),1)]/b1b2 ; % len in b1b2 coordinates
len_element = sqrt(sum(len_element.^2,2)) ;
totk = max(size(rearrklist)); PQ = zeros(totk,totk,2); MQ = PQ;

for xy = 1:2
    for kp = 1:totk
        MQ(:,kp,xy) =   rearrklist(:,xy) - rearrklist(kp,xy);
        PQ(:,kp,xy) = +(rearrklist(:,xy) + rearrklist(kp,xy));
    end
end

save([fdir,fout],'rearrbdin','rearrklist','len_element','PQ','MQ','-mat');

end