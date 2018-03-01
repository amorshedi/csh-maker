%does not support the case on one carbon


function csh_maker

inp.Ca_Si_target = 1.5;
inp.nx = 1;inp.ny = 1;
inp.finite = 0;  %1 means finite sized

nch2 = 1; %number of carbons in the polymer

force_field = 0; %1:interface 2:csh-ff 3:clay-ff

inp.wat_dist = -1.7; %distance of water from the surface
lay_gap = 0; %gap between layers


%*****************  end of input data *****************************

[lay1,~] = create_layer(1,inp); %first argument: which layer
[lay2,get_back] = create_layer(2,inp);

trans_mat = lay1.trans; trans_mat(:,3) = 2*trans_mat(:,3);
aa = lay1.trans(1,1:2);
bb = lay1.trans(2,1:2);
cc = lay1.trans(3,:)*2;

nlay1 = length(lay2.elems);
elems = [lay1.elems;lay2.elems];
coords = [lay1.coords;lay2.coords + [-2.0383    2.6077 0]]; %the difference between the coordinates of two od is added
coords = coords - get_back;
% cc = cc  + cc/norm(cc)*lay_gap;
cc(3) = cc(3)  + lay_gap;

trans_mat(:,3) = cc;


%% add polymer atoms / bonding information
if nch2>0
    [coords,elems,bonds,angles,dihedrals,trans_mat] = add_polymer2(coords,elems,nch2,trans_mat,nlay1);
else
    elems(strcmp(elems,'Os')) = {'Oh'};
    
    b_info = {{'Oc' 'Si' 2.2  1}
        {'Oh'  'Sib'  2.2  2}
        {'Oh'  'H'   1.3  3}
        {'Sib'  'Od'  2.2  4}
        {'Si'  'Ob'  2.2  5}
        {'Sib' 'Ob'  2.2  6}
        {'Ow'  'Hw'  1.3  7}
        {'Oh'  'Si'  2.2  2}
        {'Sib' 'Oh'  'H'   2.2 1}
        {'Si' 'Oh'  'H'   2.2 1}
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
    
    bonds = [(1:size(bonds,1))' bonds];
    angles = [(1:size(angles,1))' angles];
    dihedrals = [(1:size(dihedrals,1))' dihedrals];
end
% elems(strcmp(elems,'On')) = {'Od'};



% [coords,elems] = change_si_oh2(coords,elems,trans_mat);


%% Charge manipulation

[vec_ch,atom_count_id] = charge(elems,force_field); 

is_close(coords,elems,.8)


%% Lammps data file
for i = 1:length(vec_ch)
    Atom_type(atom_count_id(:,i)) = i;
    Charge_atom(atom_count_id(:,i)) = vec_ch(i);
end

[masses,pair_coeffs,bond_coeffs,angle_coeffs,dihedral_coeffs] = datafile_stuff;

molecules = zeros(size(elems));
molecules(ismember(elems,{'C','Hc'})) = 1;

aa = [aa 0]; bb = [bb 0];
[R,~] = AxelRot(subspace(aa',[1;0;0])*180/pi,[0 0 1]);
coords = coords*R; 
aa = aa*R;
bb = bb*R;
cc = cc*R;
coords = [(1:size(coords,1))' molecules Atom_type' Charge_atom' coords];

nel = length(elems);


aa = trans_mat(1,:); bb = trans_mat(2,:); cc = trans_mat(3,:);
u1 = aa/norm(aa);
u2 = bb/norm(bb);
u3 = cc/norm(cc);

lx=norm(aa);
xy=norm(bb)*dot(u1,u2); %bcos(gam) where thet is the angle between u1 and u2
xz=norm(cc)*dot(u1,u3); %u1*e3
ly=sqrt(norm(bb)^2-xy^2);
yz=(norm(bb)*norm(cc)*dot(u2,u3)-xy*xz)/(ly);
lz=sqrt(norm(cc)^2-xz^2-yz^2);

xlo = 0; xhi = lx; ylo = 0;yhi = ly; zlo = 0;zhi = lz;

% minc = min(coords(:,end-2:end));maxc = max(coords(:,end-2:end));
% xlo = minc(1)-2;  ylo = minc(2)-2;  zlo = minc(3)-2;
% xhi = maxc(1)+2;  yhi = maxc(2)+2;  zhi = maxc(3)+2;

fid = fopen('data.dat','w');
out = ['\n' num2str(nel) ' atoms\n'...
    num2str(size(bonds,1)) ' bonds\n'...
    num2str(size(angles,1)) ' angles\n'...
    num2str(size(dihedrals,1)) ' dihedrals\n\n'...
    ...
    num2str(size(masses,1)) ' atom types\n'...
    num2str(size(bond_coeffs,1)) ' bond types\n'...
    num2str(size(angle_coeffs,1)) ' angle types\n'...
    num2str(size(dihedral_coeffs,1)) ' dihedral types\n\n'...
    ...
    num2str([xlo xhi]) ' xlo xhi\n'...
    num2str([ylo yhi]) ' ylo yhi\n'...
    num2str([zlo zhi]) ' zlo zhi\n'...
    num2str([xy,xz,yz]) ' xy xz yz\n\n'...
    ...
    'Masses\n\n' outputify(masses) '\n\n'...
    ...
    'Pair Coeffs\n\n'...
    outputify(pair_coeffs) '\n\n'...
    'Bond Coeffs\n\n'...
    outputify(bond_coeffs) '\n\n'...
    'Angle Coeffs\n\n'...
    outputify(angle_coeffs) '\n\n'...
    'Dihedral Coeffs\n\n'...
    outputify(dihedral_coeffs) '\n\n'...
    ...
    'Atoms\n\n'...
    outputify(coords) '\n'...
    ...
    'Bonds\n\n'...
    outputify(bonds) '\n'...
    ...
    'Angles\n\n'...
    outputify(angles) '\n\n'...
    ...
    'Dihedrals\n\n'...
    outputify(dihedrals) '\n\n'...
    ];

fprintf(fid,out);


end


% functions *********************************************************
function [out,get_back] = create_layer(lay_num,inp)

wat_dist = inp.wat_dist;

[coords, elems, trans_mat,get_back] = atom_types (inp.nx, inp.ny,lay_num);

%Increasing Ca/Si Ratio of the TB11 layer to CSH of Ca/Si > 1
[coords,elems,num_exist] = inc_cs(coords,elems,trans_mat,inp.Ca_Si_target,lay_num);

if inp.finite==1   %% finite size
    
    % make all the silicons tetrahedral
    [coords,elems,num_exist] = make_tetra(coords,elems,trans_mat,num_exist);
    
    % remove edge bridging tetrahedra, monomers, pairing with one silicon
    [coords,elems,num_exist] = rmv_edge_tetra (coords,elems,num_exist);
    
    % fix calciums
    [coords,elems,num_exist] = fix_ca(coords,elems,trans_mat,num_exist);
end

%Oh-H is added to Cw
[coords,elems,num_exist] = add_cwoh(coords,elems,num_exist);

%Adding H to the ending oxygens in the system to balance the system charge
[coords,elems,num_exist] = change_si_oh(coords,elems,trans_mat,num_exist);


%% prepare output
num_exist = logical(num_exist);
out.coords = coords(num_exist,:);
out.elems = elems(num_exist);
out.trans = trans_mat;

end
