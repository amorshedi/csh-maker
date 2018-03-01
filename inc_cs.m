% by removing charge neutralized Si-O2 groups. If dimer removed, add H to ALL Oc.
% Removing SiO2 groups to comply with Drierketten rule (MCL=3N+2)

function [coords,elems,num_exist] = inc_cs(coords,elems,trans_mat,Ca_Si_target,lay_num)

sib = find(strcmp(elems,'Sib'));
od_id = find(strcmp(elems,'Od'));
od_extend = box_extend(coords(od_id,:),trans_mat);

od_sib = [];
for isi = sib'
    
    dist = pnt_dist(od_extend,coords(isi,:));
    tmp = find(dist<1.8);
    tmp = mod(tmp,length(od_id)); %revert back to the actual numbering
    tmp(tmp==0) = length(od_id);
    
    od_sib = [od_sib;od_id(tmp)'];
end

mid = mean(coords(sib,3));
sib_top = sib(coords(sib,3)>mid);
sib_bot = sib(coords(sib,3)<mid);
n_sib = length(sib);

num_exist = ones(length(elems),1);
num_si = sum(ismember(elems,{'Si','Sib'}));

%compute si that should be removed from natcomm paper
num_ca = sum(ismember(elems,{'Ca','Cw'}));
Ca_Si_ratio = num_ca / num_si;
n_sio2_rm = ceil(num_ca*(1/Ca_Si_ratio-1/Ca_Si_target)); %p.1 for deriv
si_check = n_sio2_rm;

% If more si than bridging should be removed, this gives you the number of bridging si
n_sio2_rm = (n_sio2_rm>n_sib)*n_sib + (n_sio2_rm<=n_sib)*n_sio2_rm; 

n_sib_top_rm = ceil(n_sio2_rm/2);
n_sib_bot_rm = n_sio2_rm - n_sib_top_rm;

%excluce the si in the middle

% mid = mean (coords);
% [~,id1] = min(pnt_dist(coords(sib_top,:),mid));
% chosen_sib1 = sib_top(id1);
% [~,id1] = min(pnt_dist(coords(sib_bot,:),coords(chosen_sib1,:)));
% chosen_sib2 = sib_bot(id1);
% if lay_num==1
%     sib_top = setdiff(sib_top,chosen_sib1);
% else
%     sib_bot = setdiff(sib_bot,chosen_sib2);
% end

rand_rm_top = randsample(sib_top,n_sib_top_rm);
rand_rm_bot = randsample(sib_bot,n_sib_bot_rm);
rand_rm = [rand_rm_top;rand_rm_bot];

num_exist (rand_rm) = 0;
num_exist (od_sib(ismember(sib,rand_rm),:)) = 0;
if si_check>n_sib %dimers hould be removed
    
    %list of top and bottom dimers
    sic = Si_neighb(Si_neighb(:,7)==1,1); %id of Si near Oc
    si_extend = cellfun(@(x)coords(sic,:)+x,im_cart,'uni', 0);
    si_extend = cell2mat(si_extend);
    chosen = []; %chosen is eventually dimer pairs
    for i = 1:length(sic)
        dist = sqrt(sum((si_extend-coords(sic(i),:)).^2,2));
        
        tmp = mod(find(dist<3.5),length(sic))';
        tmp(tmp==0) = length(sic);
        chosen = [chosen;sic(tmp)'];
        
    end
    
    %remove 1-2 and 2-1 cases
    [~,b] = unique(sort(chosen,2),'rows');
    chosen = chosen(b,:);
    
    mid = mean(coords(chosen(:,1),3));
    idx = coords(chosen(:,1),3)>mid; %choose the top and bot dimers
    dim_top = chosen(idx,:);
    dim_bot = chosen(~idx,:);
    
    n_rm = ceil((si_check - n_sib)/2); %just count num of top and bot dimers that should be removed
    n_rm_top = ceil(n_rm/2);
    n_rm_bot = n_rm - n_rm_top;
    
    rand_rm_top = dim_top(randsample(size(dim_top,1),n_rm_top),:); %choose some dimers randomly
    rand_rm_bot = dim_bot(randsample(size(dim_bot,1),n_rm_bot),:);
    rand_rm = [rand_rm_top;rand_rm_bot];
    
    num_exist (rand_rm(:)) = 0; 
    o_id = Si_neighb(ismember(Si_neighb(:,1),rand_rm(:)),3:6); %remove Obs of the dimer
    idx = strcmp(elems(o_id(:)),'Ob');
    num_exist (o_id(idx)) = 0;
    
    elems(o_id(~idx)) = {'Oh'}; %saturate the Ocs of the dimer
    elems(end+1:end+sum(~idx)) = {'H'};
    num_exist (end+1:end+sum(~idx)) = 1;
    coords(end+1:end+sum(~idx),:) = coords(o_id(~idx),:) + [0 1 0];
end