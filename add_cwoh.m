%randomly choose cw, add oh - h to them

function [coords,elems,num_exist] = add_cwoh(coords,elems,num_exist)

nelem_layer = elems(logical(num_exist)); %only the elements in the current model
num_si = sum(ismember(nelem_layer,{'Si','Sib'}));
num_ca = sum(ismember(nelem_layer,{'Ca','Cw'}));
num_H = sum(strcmp(nelem_layer,'H'));
ca_si_ratio = num_ca / num_si;

num_add_cw_OH = round(num_si * (ca_si_ratio - 1)-num_H); %-num_H because we are counting the OcH as CaOH
num_add_OH_top = floor(num_add_cw_OH/2);
num_add_OH_bot = num_add_cw_OH-num_add_OH_top;

mid = mean(coords);
cw_top = find(strcmp(elems,'Cw') & coords(:,3)>mid(3));
cw_bot = setdiff(find(strcmp(elems,'Cw')),cw_top);

samp_top = randsample(cw_top,num_add_OH_top);
samp_bot = randsample(cw_bot,num_add_OH_bot);

oh_unit = [1 1 0;-1 1 0]/2;
coords = [coords;cell2mat(cellfun(@(x)oh_unit+x,num2cell(coords([samp_top;samp_bot],:),2),'uni',0))];
elems = [elems;repmat({'Oh';'H'},num_add_OH_top+num_add_OH_bot,1)];

num_exist = [num_exist;ones(2*(num_add_OH_top+num_add_OH_bot),1)];