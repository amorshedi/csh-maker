function [coords,elems,num_exist] = make_tetra(coords,elems,trans_mat,num_exist)
%if you consider images and look at the silicons in the main box, you can
%make all Si tetrahedral. Some of these oxygen will be from the images

si_id = find(ismember(elems,{'Si','Sib'})); %id of silicons

o_id = find(ismember(elems,{'Oc','Od','Ob'}));
o_extend = box_extend(coords(o_id,:),trans_mat); %O=coordinates of diff types of O atoms

rcut = 2.;
image_oxygens = [];
for i = si_id'
    dist = pnt_dist(o_extend,coords(i,:));
    tmp = find(dist<rcut);    
    
%a few of the oxygens around the current si may be in the image
    image_oxygens = [image_oxygens;setdiff(tmp,tmp((tmp>13*length(o_id))&(tmp<=14*length(o_id))))];

end
chosen = mod(image_oxygens,length(o_id)); %revert back to the actual numbering
chosen(chosen==0) = length(o_id);
o_id = o_id(chosen);

%add the images
coords = [coords;o_extend(image_oxygens,:)];
elems = [elems;elems(o_id)];
num_exist = [num_exist;ones(length(o_id),1)];

%remove oc and od (because Ob makes two Si tetrahedral)
idx = ismember(elems(o_id),{'Oc','Od'});
num_exist(o_id(idx)) = 0;