function [coords,elems,num_exist] = rmv_edge_tetra (coords,elems,num_exist)

% determine the Ob connected to only one Si
ob_id = find(ismember(elems,'Ob'));
si_id = find(ismember(elems,{'Si','Sib'})); %id of silicons
si_coords = coords(si_id,:); 

rcut = 1.8; end_si = []; ob_one_si = [];
for iob = ob_id'
    dist = pnt_dist(si_coords,coords(iob,:));
    tmp = find(dist<rcut); 
    
    if length(tmp)==1
        end_si = [end_si;si_id(tmp)];
        ob_one_si = [ob_one_si;iob];
    end
end

%remove the ending bridging tetrahedra at the end of chains
%remove monomers
o_id = find(ismember(elems,{'Oc','Od','Ob'}));
o_coords = coords(o_id,:); %O=coordinates of diff types of O atoms
rcut = 1.8; rm_list = [];
for i = end_si'
    dist = sqrt(sum((o_coords-coords(i,:)).^2,2));
    tmp = o_id(dist<rcut); %id of the oxygens of the tetrahedra (always 4 after modifications done earlier)
    
    next_ob = tmp(strcmp(elems(tmp),'Ob')&~ismember(tmp,ob_one_si)); %the Ob not on the edge
    if isempty(next_ob)
        rm_list = [rm_list;[tmp;i]];%the Ob not on the edge
        continue
    end
    
    if ismember('Od',elems(tmp)) %if the first Si is bridging
        rm_list = [rm_list;[i;setdiff(tmp,next_ob)]]; %I want to remove 2*Od, 1*Si and the Ob on the edge
        
    else %If the next tetrahedra is dangling, remove both the edge and the dangling
        
        dist = pnt_dist(si_coords,coords(next_ob,:));
        next_si = setdiff(si_id(dist<rcut),i);
        dist = pnt_dist(o_coords,coords(next_si,:));
        tmp2 = o_id(dist<rcut); %id of the oxygens of the tetrahedra
        
        if ismember('Od',elems(tmp2)) %if the next si after side is bridging
            rm_list = [rm_list;[tmp;i]]; %remove all of first tetrahedron
            rm_list = [rm_list;[tmp2(strcmp(elems(tmp2),'Od'));next_si]]; %remove the remaining of the second except Ob
        end
        
    end
    
end

rm_list = unique(rm_list); %I think a chain of three Si causes duplicates. No big deal
num_exist(rm_list) = 0;


