function [xlo, klist, full_qty, hex_vert] = reproduce_hexagon(b1b2klist,b1b2_qty,b1b2, hexagon_vert)

b1b2klist = b1b2klist * b1b2;
hex_vert = hexagon_vert * b1b2 ; ones_vert = ones(max(size(hex_vert)),1) ;
xlo = min(hex_vert) ; xlo = xlo(1)*1.1;

in = inpolygon(b1b2klist(:,1),b1b2klist(:,2),hex_vert(:,1),hex_vert(:,2));
klist = b1b2klist(in,:); full_qty = reshape(b1b2_qty(in),[],1);

%translate from up to down by b2 vector
disp_b2 = b1b2(2,:) ; % this is the b2 vector along Y-axis
new_hex_vert = hex_vert + [ disp_b2(1)*ones_vert, disp_b2(2)*ones_vert];
in = inpolygon(b1b2klist(:,1),b1b2klist(:,2),new_hex_vert(:,1),new_hex_vert(:,2));
ones_in = ones(max(size( b1b2klist(in,:) )),1) ;
klist = [klist; b1b2klist(in,:)- [ disp_b2(1)*ones_in, disp_b2(2)*ones_in] ];
full_qty = [full_qty; reshape(b1b2_qty(in),[],1) ];

%translate from right bottom to left up corner by b1
disp_b1 = b1b2(1,:); % this is the b1 vector 
new_hex_vert = hex_vert + [ disp_b1(1)*ones_vert, disp_b1(2)*ones_vert];
in = inpolygon(b1b2klist(:,1),b1b2klist(:,2),new_hex_vert(:,1),new_hex_vert(:,2));
ones_in = ones(max(size( b1b2klist(in,:) )),1) ;
klist = [klist; b1b2klist(in,:)- [ disp_b1(1)*ones_in, disp_b1(2)*ones_in] ];
full_qty = [full_qty; reshape(b1b2_qty(in),[],1) ];

%translate from up to down by b1+b2 vector
disp_b1b2 = b1b2(2,:)+b1b2(1,:) ; % this is the b1+b2 vector along Y-axis
new_hex_vert = hex_vert + [ disp_b1b2(1)*ones_vert, disp_b1b2(2)*ones_vert];
in = inpolygon(b1b2klist(:,1),b1b2klist(:,2),new_hex_vert(:,1),new_hex_vert(:,2));
ones_in = ones(max(size( b1b2klist(in,:) )),1) ;
klist = [klist; b1b2klist(in,:)- [ disp_b1b2(1)*ones_in, disp_b1b2(2)*ones_in] ];
full_qty = [full_qty; reshape(b1b2_qty(in),[],1) ];

end