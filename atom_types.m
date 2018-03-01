%#
function  [final_coords, elems, trans_mat,get_back] = atom_types (nx, ny, lay_num)
% This function finds the atom types in the structure of calcium-alumino-silcates
% and constructs a supercell plate of Tobermorite

% nx = 2; ny = 1;lay_num=2;

[elems, coords, trans_mat] = coordinate_transformation;

%remove water
idx = ismember(elems,{'Ow','Hw'});
coords(idx,:) = [];
elems(idx)=[];

aa=trans_mat(1,:);
bb=trans_mat(2,:);
cc=trans_mat(3,:);

%distinguishing Cw
cw_top = strcmp(elems,'Ca')&(coords(:,3)<19.73 & coords(:,3)>14.89);
cw_bot = strcmp(elems,'Ca')&(coords(:,3)<7.83 & coords(:,3)>2.67);
idx = cw_top|cw_bot;
elems(idx) = {'Cw'};
cw_bot = find(cw_bot);
cw_top = find(cw_top);

%select the layer in the middle
idx_bot = find((~strcmp(elems,'Cw'))&(coords(:,3)<17.5 & coords(:,3)>5.));
idx_bot = [idx_bot;[cw_top(1:2:end);cw_bot(1:2:end)]];
elems1 = elems(idx_bot);
coords1 = coords(idx_bot,:);

%the other layer
idx_top = setdiff(1:length(elems),idx_bot)';
elems2 = elems(idx_top);
coords2 = coords(idx_top,:);

% seperating the layers (see p. 5 of the notebook)
mid = mean(coords(:,3));
bot_bot = find(coords2(:,3)<mid); %some of the layer is at the bottom and some is at top of the box
coords2(bot_bot,:) = coords2(bot_bot,:) + cc; %translating in the c direction

get_back = cc*5/cc(3); %how much empty space there is between bottom of abc triad and material

if lay_num == 1
    elems = elems1;
    coords = coords1;
else
    elems = elems2;
    coords = coords2;
end

trans_mat(3,:) = trans_mat(3,:)/2;

%oxygens: Od or Oc?
extend = box_extend(coords,trans_mat);
o_id = find(strcmp(elems,'O'));
rcut = 2.8;
for io = o_id' %if there is a Ca around O -> Oc and if Cw -> Od
    dist = pnt_dist(extend,coords(io,:));
    tmp = find(dist<rcut);
    tmp = mod(tmp,size(coords,1)); %revert back to the actual numbering
    tmp(tmp==0) = size(coords,1);

    if ismember('Ca',elems(tmp))
        elems{io} = 'Oc';
    elseif sum(ismember(elems(tmp),'Si'))>1
        elems{io} = 'Ob';
    end
end
elems(strcmp(elems,'O')) = {'Od'};

%differentiate bridging Si
si_id = find(strcmp(elems,'Si'));
od_id = strcmp(elems,'Od');
od_extend = box_extend(coords(od_id,:),trans_mat);

for isi = si_id'
    
    dist = pnt_dist(od_extend,coords(isi,:));
    
    if sum(dist<2)>1
        elems(isi) = {'Sib'};
    end
     
end


% constructing a n_x * n_y supercell of Tobermorite Layer
[im_a,im_b] = meshgrid(0:nx-1,0:ny-1);
im_cart = [im_a(:),im_b(:)]*trans_mat(1:2,1:2); 
im_cart(:,3) = 0;

dim1 = ones(1,size(im_cart,1));
im_cart = mat2cell(im_cart,dim1,3);

final_coords = cellfun(@(x)coords+x,im_cart,'uni',0);
final_coords = cell2mat(final_coords);

elems = repmat(elems,nx*ny,1);

trans_mat(1,:) = trans_mat(1,:) * nx;
trans_mat(2,:) = trans_mat(2,:) * ny;
trans_mat(3,:) = trans_mat(3,:); %times 2 needed if writing for two layers 

end