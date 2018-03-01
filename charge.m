function [vec_ch,atom_count_id] = charge(elems,force_field)

%             1     2    3    4   5     6   7     8    9   10   11   12  13
all_atoms = {'Oc' 'Ca' 'Oh' 'Od' 'Ob' 'Cw' 'Si' 'Sib' 'Ow' 'Hw' 'H' 'Hc' 'C'};
atom_count_id = cell2mat(cellfun(@(x)strcmp(elems,x),all_atoms,'uni',0));
atom_count = sum(atom_count_id);

vec_ch = [-0.94,1.500,-0.6700,-1.04,-0.58,1.7,1 ,1,-0.82,0.41,0.4, 0.1, -0.2]; %INTERFACE-CVFF
total_charge = sum(atom_count.*vec_ch);

all_o_share=total_charge/sum(atom_count([1 3 4 5]));
vec_ch([1 3 4 5]) = vec_ch([1 3 4 5]) - all_o_share;


check_charge = sum(atom_count.*vec_ch); %last term was calculated