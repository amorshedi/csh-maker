% difference with add_polymer1: this adjusts distance automatically
% choose two close od instead of silicon
% exclude the other Od of the silicon that is chosen

function [coords,elems,bonds,angles,dihedrals,trans_mat] = add_polymer2(coords,elems,nch2,trans_mat,nlay1)

sio3_0 = [1.076 -.517 -3.726; 1.977 .799 -4.152; 2.016 -1.821 -4.099; -.192 -.61 -4.78];
sio3_1 = [.799 2.544 5.31; -.272 1.354 5.687; 2.033 2.309 6.382; 0.01447 3.928 5.7578];

ch2_1 = [.545 -.528 -1.982; -.555 -.449 -1.960; .776 -1.525 -1.562];
ch2_2 = [1.167 .559 -1.108; .935 1.557 -1.521; 2.268 .477 -1.112];


ch2_trans = [.649 .473 .353] - ch2_1(1,:); %first one is the coords of the third C
si_trans = ch2_2(1,:) - sio3_0(1,:)-ch2_trans/2;

pln = createPlane(sio3_0(1,:), ch2_trans/norm(ch2_trans));

sio3_2 = projPointOnPlane(sio3_0,pln);
sio3_2 = sio3_2 - (sio3_0-sio3_2);
% sio3_3 = sio3_1 - (sio3_1(1,:)-sio3_0(1,:));

pol_elems = repmat({'C';'Hc';'Hc'},nch2,1);
pol_elems = [{'Sib';'Ob';'Ob';'Od'};pol_elems;{'Sib';'Ob';'Ob';'Od'}];

pol_coords = [];
for ich2 = 1:nch2

    pol_coords = [pol_coords;ch2_trans*(ceil(ich2/2)-1)+mod(ich2,2)*ch2_1 + ~mod(ich2,2)*ch2_2];
    
end

pol_coords = [sio3_0;pol_coords];

if mod(nch2,2)==1
    pol_coords = [pol_coords;sio3_2+2*si_trans+ch2_trans*(ceil(ich2/2)-1)];
else
    pol_coords = [pol_coords;sio3_1+ch2_trans*(ceil(ich2/2)-3)];
end

b_info = {
    {'Sib' 'C'   2    8}
    {'C'   'Hc'  2    9}
    {'C'   'C'   2    10}
    {'Od'   'Sib' 'C'   2.5  13}
    {'Ob'   'Sib' 'C'   2.5  13}
    {'Hc'   'C'  'Hc'   2.  14}
    {'C'    'C'  'Hc'   2.  15}
    {'Sib'   'C' 'Hc'   2.  16}
    {'Sib'   'C' 'C'   2.  17}
    {'C'   'C' 'C'     2.  18}
    {'Od'     'Sib'   'C'    'Hc'  2  1}
    {'Od'     'Sib'   'C'    'C'   2  2}
    {'Ob'     'Sib'   'C'    'Hc'  2  3}
    {'Ob'     'Sib'   'C'    'C'   2  4}
    {'Hc'     'C'     'C'    'Hc'  2  5} 
    {'C'     'C'     'C'    'Hc'  2  6} 
    {'Sib'   'C'     'C'    'Hc'  2  7} 
    {'Sib'   'C'     'C'    'C'   2  8}
    {'C'     'C'     'C'    'C'   2  9}
    };

trans_mat1 = [30 0 0;0 30 0;0 0 30];

[pol_bonds,pol_angles,pol_dihedrals] = form_bonds(pol_coords,pol_elems,trans_mat1,b_info);
% 
% elems = pol_elems;
% coords = pol_coords;
% vmd

pol_elems = pol_elems(5:end-4);
pol_coords = pol_coords(5:end-4,:);

%adjust interlayer distance
lay_shift = (-0.9272+norm(ch2_trans)/2*nch2)/trans_mat(3,3)*trans_mat(3,:);
coords(nlay1+1:end,:) = coords(nlay1+1:end,:) + lay_shift;
trans_mat(3,:) = trans_mat(3,:) + 2*lay_shift;


%choose two Si atoms
mid = mean(coords);
top_sib_id = find(strcmp(elems,'Sib') & coords(:,3)>mid(3));
bot_sib_id = find(strcmp(elems,'Sib') & coords(:,3)<mid(3));

[~,id1] = min(pnt_dist(coords(bot_sib_id,:),mid));
chosen_sib1 = bot_sib_id(id1); 
[~,id1] = min(pnt_dist(coords(top_sib_id,:),coords(chosen_sib1,:)));
chosen_sib2 = top_sib_id(id1); 


% %find the Od and Ob of each the Si
rcut = 2;
od1_id = find(pnt_dist(coords,coords(chosen_sib1,:))<rcut & ismember(elems,{'Od','Os'}));
od2_id = find(pnt_dist(coords,coords(chosen_sib2,:))<rcut & ismember(elems,{'Od','Os'}));

cmod1 = coords(od1_id(1),:);
cmod2 = coords(od2_id(1),:);

%coordinate of the end C of polymer
cpc1 = pol_coords(1,:);

%move the first C of polymer to model's first chosen Od
pol_coords = pol_coords + (cmod1 - cpc1);

%update
cpc1 = pol_coords(1,:); 
cpc2 = pol_coords(end-2,:);

mod_vec = cmod2 - cmod1; %vector connecting model's od
ch2_trans = cpc2 - cpc1; %polymer vector

%stretch polymer
fact1 = norm(mod_vec)/norm(ch2_trans) - 1; 
pol_plane = createPlane(cmod1,ch2_trans/norm(ch2_trans));
for i = 1:size(pol_coords,1)
    ax_vec = distancePointPlane(pol_coords(i,:),pol_plane)*ch2_trans/norm(ch2_trans)*fact1;
    pol_coords(i,:) = pol_coords(i,:) + ax_vec;
end


ang = dot(mod_vec,ch2_trans)/(norm(ch2_trans)*norm(mod_vec));
rot_ax = cross(mod_vec,ch2_trans);

[rot,~] = AxelRot(-acosd(ang),rot_ax/norm(rot_ax),cmod1); %I think u does not have to be a unit vector

pol_si_vec = pol_coords - cmod1;

%transformed polymer
pol_coords = (rot*pol_si_vec')' + cmod1;

coords = [coords;pol_coords];
elems = [elems;pol_elems];

% if the oxygen replaced with C is an Os (Oh connected to Si), remove its H
rm_list = [od1_id(1) od2_id(1)];
for i = rm_list
    if strcmp(elems(i),'Os')
        tmp = pnt_dist(coords,coords(i,:))<1.1;
        rm_list = [rm_list find(strcmp(elems,'H') & tmp)];
    end
end

%because I don't have Oh-Si-C angle coeffs at the moment, if the other
%oxygen is Os, change it to Od and remove its oxygen.
os_list = [od1_id(2) od2_id(2)];
for i = os_list
    if strcmp(elems(i),'Os')
        tmp = pnt_dist(coords,coords(i,:))<1.1;
        rm_list = [rm_list find(strcmp(elems,'H') & tmp)];
    end
end
elems(os_list) = {'Od'};
    
coords(rm_list,:) = [];
elems(rm_list) = [];
elems(strcmp(elems,'Os')) = {'Oh'};


b_info = {{'Oc' 'Si' 2.2  1}
    {'Oh'  'Sib'  2.2  2}    
    {'Oh'  'H'   1.3  3}
    {'Sib'  'Od'  2.2  4}
    {'Si'  'Ob'  2.2  5}
    {'Sib' 'Ob'  2.2  6}
    {'Ow'  'Hw'  1.3  7}
    {'Oh'  'Si'  2.  2}
    {'Sib' 'Oh'  'H'   2.2 1}
    {'Si'  'Ob'  'Si'  2.2 2}
    {'Si'  'Ob'  'Sib' 2.2 3}
    {'Oc'   'Si'  'Oc'  2.2  4}
    {'Oc'   'Si'  'Ob'  2.2  5}
    {'Ob'   'Si'  'Ob'  2.2  6}
    {'Oh'   'Sib' 'Oh'  2.2  7}
    {'Oh'   'Sib' 'Ob'  2.2  8}
    {'Ob'   'Sib' 'Ob'  2.2  9}
    {'Hw'   'Ow'  'Hw'  1.3  10}
    {'Oh'   'Sib' 'Od'  2.2  11}
    {'Od'   'Sib' 'Ob'  2.2  12}
    };

[bonds,angles,dihedrals] = form_bonds(coords,elems,trans_mat,b_info);


%choose two Si atoms
% mid = mean(coords);
top_sib_id = find(strcmp(elems,'Sib') & coords(:,3)>mid(3));
bot_sib_id = find(strcmp(elems,'Sib') & coords(:,3)<mid(3));

[~,id1] = min(pnt_dist(coords(bot_sib_id,:),mid));
chosen_sib1 = bot_sib_id(id1); 
[~,id1] = min(pnt_dist(coords(top_sib_id,:),coords(chosen_sib1,:)));
chosen_sib2 = top_sib_id(id1); 


% %find the Od and Ob of each the Si
rcut = 2;
od1_id = find(pnt_dist(coords,coords(chosen_sib1,:))<rcut & strcmp(elems,'Od'));
od2_id = find(pnt_dist(coords,coords(chosen_sib2,:))<rcut & strcmp(elems,'Od'));
ob1_id = find(pnt_dist(coords,coords(chosen_sib1,:))<rcut & strcmp(elems,'Ob'));
ob2_id = find(pnt_dist(coords,coords(chosen_sib2,:))<rcut & strcmp(elems,'Ob'));

current_id = 1:length(pol_elems) + 8; %8 is cause of two silicon heads
model_id = [chosen_sib1 ob1_id' od1_id' length(elems)+(1-length(pol_elems):0) chosen_sib2 ob2_id' od2_id];

btmp = pol_bonds(:,2:end);
btmp1 = btmp;
atmp = pol_angles(:,2:end);
atmp1 = atmp;
dtmp = pol_dihedrals(:,2:end);
dtmp1 = dtmp;
for id = current_id
    
    btmp1(btmp==id) = model_id(id);
    atmp1(atmp==id) = model_id(id);
    dtmp1(dtmp==id) = model_id(id);
    
end

btmp1 = [pol_bonds(:,1) btmp1];
atmp1 = [pol_angles(:,1) atmp1];
dtmp1 = [pol_dihedrals(:,1) dtmp1];

bonds = [bonds; btmp1];
angles = [angles; atmp1];
dihedrals = [dihedrals; dtmp1];

bonds = [(1:size(bonds,1))' bonds];
angles = [(1:size(angles,1))' angles];
dihedrals = [(1:size(dihedrals,1))' dihedrals];












