%move the calcium that is unhappy in the nonperiodic to its image
%so that it can be happy there with Oc
function [coords,elems,num_exist] = fix_ca(coords,elems,trans_mat,num_exist)


ca_id = find(strcmp(elems,'Ca') & num_exist); 
ca_coords = coords(ca_id,:); 

oc_id = find(ismember(elems,'Oc') & num_exist);
oc_coords = coords(oc_id,:); 

rcut = 2.7; ca_bad = [];
for ica = ca_id'
    dist = pnt_dist(oc_coords,coords(ica,:));
    tmp = find(dist<rcut);
    
    if length(tmp)==1
        ca_bad = [ca_bad;ica]; %calciums that are undercoordinated
    end
end

ca_extend = box_extend(ca_coords,trans_mat);

rcut = 2.7;
image_ca = []; sz = length(ca_id); box_ca = [];
for ioc = oc_id'
    dist = pnt_dist(ca_extend,coords(ioc,:));
    tmp = find(dist<rcut); 
    
    tmp2 = tmp(tmp<=13*sz | tmp>14*sz);
    tmp = mod(tmp2,length(ca_id));
    tmp(tmp==0) = length(ca_id);
    tmp = ca_id(tmp); %id of calcium whose image make an oc happy
    
    image_ca = [image_ca;tmp2(ismember(tmp,ca_bad))];
    box_ca = [box_ca;tmp(ismember(tmp,ca_bad))];

end
[box_ca,idx] = unique(box_ca);
image_ca = image_ca(idx);

%add the ca from image
coords = [coords;ca_extend(image_ca,:)];
elems = [elems;elems(box_ca)];
num_exist = [num_exist;ones(length(box_ca),1)];

num_exist(box_ca) = 0;





