%change a list of oxygens that are connected to only one silicon to Oh and
%add an H

function [coords,elems] = change_si_oh2(coords,elems,trans_mat)

%find o connected to one Si
o_id = find(ismember(elems,{'Ob','Od'}));
si_id = find(ismember(elems,{'Si','Sib'}));
si_extend = box_extend(coords(si_id,:),trans_mat);
rcut = 2; o_one_si = [];
for io = o_id'
    
    dist = pnt_dist(si_extend,coords(io,:));
    
    if sum(dist<rcut)==1
        o_one_si = [o_one_si; io];
    end
    
end

if ~isempty (o_one_si)
    num_si = length(si_id);
    num_ca = sum(ismember(elems,{'Ca','Cw'}));
    ca_si_ratio = num_ca / num_si;
    
    mid = mean(coords);
    o_one_top = o_one_si(coords(o_one_si,3)>mid(3));
    o_one_bot = o_one_si(coords(o_one_si,3)<mid(3));
    
    n_o_one = length(o_one_si);
    n_sioh_add = round(num_si*(ca_si_ratio-.46));
    if n_sioh_add>n_o_one
        o_rand_top=o_one_top;
        o_rand_bot=o_one_bot;
    else
        n_sioh_top=round(n_sioh_add/2);
        n_sioh_bot=n_sioh_add-n_sioh_top;
        o_rand_top=randsample(o_one_top,n_sioh_top); %some dangling Ob
        o_rand_bot=randsample(o_one_bot,n_sioh_bot);
    end
    
    o_change = [o_rand_top;o_rand_bot];
    
    elems(o_change)={'Oh'};
    elems = [elems;repmat({'H'},length(o_change),1)];
    coords = [coords;coords(o_change,:)+[0 1 0]];
end