%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kz= 0 fermi sheet contour plot to find scattering vectors between fermi sheets
% Shinibali; July 20,2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [frac_kf,full_bandind,len_dkf] = kpcontour(orbitals,Ef,totkf_in_ring,b1b2,hexagon_vert)

nx = 99; bzedge = 1; kx = (0 : bzedge/nx : bzedge); kx=kx(1:end-1) ;
[KX,KY]=meshgrid(kx) ; [nwan,Ek] = eigenThomale(KX(:),KY(:),orbitals);
Ek = reshape(Ek,nwan,nx,nx); Ek=Ek-Ef;

A = cell(nwan,1);

figure
hold on; box on;
for i = 1:nwan
    cont(:,:) = Ek(i,:,:);
    [A{i},~] = contour(KX,KY,cont,[0 0]);
end
hold off; set(gcf,'Visible','off');
counter = 0;
for n=1:nwan
    if isempty(A{n}) ~= 1
        totks = max(size(A{n})); counter = counter+1;
        
        rep = 1;
        nks(rep) = A{n}(2,1);
        
        while (nks(rep) ~= 0)
            colstart = 0; colend = 0;
            for s = 1:rep
                if (s ~= 1)
                    colstart = nks(s-1) + colstart;
                end
                colend = nks(s) + colend;
            end
            
            colstart = rep + 1 + colstart;
            colend = rep + colend;
            
            if (rep == 1)
                finale{rep,n} = A{n}( : , colstart : colend );
            else
                finale{rep,n} = horzcat(finale{rep-1,n}, A{n}( : , colstart : colend ));
            end
            
            if (colend == totks)
                break
            end
            
            rep = rep+1;
            nks(rep) = A{n}(2,colend+1);
            
        end
        
        finale{rep,n} = transpose(finale{rep,n});
        kband{n} = [ finale{rep,n}(:,1), finale{rep,n}(:,2)];
        totband(n) = max(size(kband{n}));
        bandij{n} = n .* ones(totband(n),1) ;
        
        if ( counter == 1)
            klist = kband{n}; bandind = bandij{n};
        else
            klist = vertcat(klist,kband{n});
            bandind = vertcat(bandind,bandij{n});
        end
        
    end
end

%% here starts the manipulation of unevenly distr kf points:

hex_angle = 2*pi/6 ;
[~, klist, bandind ] = ...
    reproduce_hexagon(klist, bandind ,b1b2, hexagon_vert) ; % reproducing the hexagonal BZ distr
theta_kf = atan2(klist(:,2),klist(:,1));  % finding angle theta for kx,ky

lt_set = theta_kf <= hex_angle; gt_set = theta_kf >= 0 ; % only selecting points within 0 and pi/3
set_kf = klist( lt_set == gt_set, :) ; set_r = sqrt(sum(set_kf.^2,2));
set_theta = atan2(set_kf(:,2),set_kf(:,1)); set_band = bandind( lt_set == gt_set );

[~,ind] = sort( set_theta ); % sorting theta
set_r = set_r(ind) ; set_kf = set_kf(ind,:) ; set_band = set_band(ind) ;

mid_kf = ( set_kf(1:end-1,:) + set_kf(2:end,:) )/2 ; % finding mid points betwn the two rings
mid_r = sqrt(sum(mid_kf.^2,2)); mid_theta = atan2(mid_kf(:,2),mid_kf(:,1));

% changing the polymonial order below can affect the fitting
p = polyfit(mid_theta,mid_r,7); % fitting a poly line of order 7, to divide the two rings
mid_theta = (0: 2*hex_angle/max(size(mid_kf)) : hex_angle)'; mid_r = polyval(p,mid_theta);
mid_kf = [ mid_r.*cos(mid_theta), mid_r.*sin(mid_theta)] ; % finding kx,ky for the separator line

fwd_theta = mid_kf(2:end,:) - mid_kf(1:end-1,:);
fwd_theta = atan2(fwd_theta(:,2),fwd_theta(:,1)); fwd_theta = [fwd_theta ; fwd_theta(end) ];
bck_theta = mid_kf(1:end-1,:) - mid_kf(2:end,:);
bck_theta = atan2(bck_theta(:,2),bck_theta(:,1)); bck_theta = [bck_theta ; bck_theta(end) ];

right_kf = set_kf(1,:); left_kf = right_kf; ones_kf = ones(max(size( mid_kf )), 1 );
right_band = set_band(1); left_band = right_band;
for ik = 1:numel(set_r)
    dist = sqrt(sum( [ [ set_kf(ik,1).*ones_kf, set_kf(ik,2).*ones_kf ] - mid_kf ].^2,2));
    [~,ind] = min(dist); % finding min dist betwn line and each kf
    dist = set_kf(ik,:) - mid_kf(ind,:) ; % finding x,y of that dist
    left_or_right = atan2( dist(2),dist(1) ); % finding corresponding theta wrt the midpoint line
    lt_set = left_or_right <= fwd_theta(ind); % whether theta < forward midpoint theta
    gt_set = left_or_right >= bck_theta(ind); % whether theta > backward midpoint theta
    
    if lt_set && gt_set
        right_kf = [right_kf; set_kf(ik,:)] ; % goes to right side of line
        right_band = [right_band; set_band(ik)] ; % goes to right side of line
    else
        left_kf = [left_kf; set_kf(ik,:)] ; % goes to left side of line
        left_band = [left_band; set_band(ik)] ; % goes to right side of line
    end
end

hex_pie_kf{1} = right_kf(2:end,:) ; % saving in cell format for the pi/3 chunk
hex_pie_kf{2} = left_kf(2:end,:) ; % excluding the false first point
right_band = right_band(2:end) ; left_band = left_band(2:end) ; % excluding the false first point

right_theta = atan2(right_kf(:,2),right_kf(:,1)); new_right_r = polyval(p,right_theta);
diff_right_kf = norm( [ new_right_r.*cos(right_theta), new_right_r.*sin(right_theta)] - right_kf );
if diff_right_kf < 0.1       % threshold to check if the line overlaps with the ring itself
    hex_pie_kf{1} = set_kf ; % this is to avoid double counting of single ring
    hex_pie_kf{2} = [] ;     % this is to avoid double counting of single ring
    right_band = set_band ;
end

frac_kf = [0,0]; len_dkf = 0; full_bandind = 0;
Rot = @(angle) [cos(angle) -sin(angle); sin(angle) cos(angle)] ; % angle in radians of pi
for iLR = 1:2     % for left or right kf
    if ~isempty(hex_pie_kf{iLR})
        xy = hex_pie_kf{iLR}; list_xy=xy;
        for irot = 1:5
            rotated_kf = Rot( hex_angle*irot ) * xy' ;
            list_xy = [list_xy; rotated_kf' ];
        end
        r_kf = sqrt(sum(list_xy.^2,2)); theta_kf = atan2(list_xy(:,2),list_xy(:,1));
        [sorted_theta,ind] = sort( theta_kf ); sorted_r = r_kf(ind) ;
        even_theta = (-pi: 2*pi/totkf_in_ring : pi)';  even_theta = even_theta(1:end-1); % evenly distr. theta
        even_r = interp1(sorted_theta,sorted_r,even_theta,'linear','extrap'); % linear interp for r_kf for given theta
        even_kf = [even_r.*cos(even_theta),even_r.*sin(even_theta) ]; % finding kx,ky from r, theta
        if iLR == 1
            full_bandind = [full_bandind; unique(right_band) * ones(size(even_kf,1),1)] ;
        else
            full_bandind = [full_bandind; unique(left_band) * ones(size(even_kf,1),1)] ;
        end
        rearr_kf = even_kf/b1b2; % finding fractional coordinates in units of b1,b2
        %shifting to square BZ
        indices  = rearr_kf(:,1) < 0 ; 
        rearr_kf( indices ,1) = rearr_kf( indices,1) + ones(max(size(rearr_kf( indices,1) ) ) ,1); 
        %shifting to square BZ
        indices  = rearr_kf(:,2) < 0 ;
        rearr_kf( indices ,2) = rearr_kf( indices,2) + ones(max(size(rearr_kf( indices,1) ) ) ,1);
        
        frac_kf = [frac_kf;rearr_kf] ;
        len_dkf = [len_dkf; even_r * ( even_theta(2)-even_theta(1) ) ] ; % length element r*dtheta
        
    end
end
frac_kf = frac_kf(2:end,:) ; len_dkf = len_dkf(2:end,:) ; full_bandind = full_bandind(2:end) ;

end