function [elems, coords, trans_mat] = coordinate_transformation
% Transforms the cartisian to fractional coordinates and vice versa and
% changing from hexagonal to orthohedral unit cell.


% Reading the text files. first line is: alen, blen, clen, alph, bet,gam
% Subsequent lines are Atom label and X, Y and Z.
fid = fopen('TB11.txt','r');
cell_dims = str2num(fgets(fid));

mydata = textscan(fid,'%s %f %f %f');
coords = cell2mat(mydata(:,2:4));
Natom = size(coords,1);
elems = mydata(:,1);
elems = elems{1};

alen  = cell_dims(1); blen = cell_dims(2); clen = cell_dims(3);
alph  = cell_dims(4) * (pi/180);
bet = cell_dims(5) * (pi/180);
gam  = cell_dims(6) * (pi/180);

% Constructing the transformation matrix in C-S-H.
trans_mat = zeros(3,3);
vol = alen * blen * clen * sqrt(1 - cos(alph)^2 - cos(bet)^2 - cos(gam)^2 + 2 * cos(alph) * cos(bet) * cos(gam));
trans_mat (1,1) = alen;
trans_mat (2,1) = blen * cos(gam);
trans_mat (2,2) = blen * sin(gam);
trans_mat (3,1) = clen * cos(bet);
trans_mat (3,2) = clen * (cos(alph) - cos(bet) * cos(gam)) / sin(gam);
trans_mat (3,3) = vol / (alen * blen * sin(gam));

% Transforming the highly tilted triclinic to nearly orthogonal unit cell
% This similar to hexagonal to orthohedral transformation technique

% Step 1: change the center of parallepiped in ab plane
coords(:,1:2) = coords(:,1:2) - 0.5;

% Step 2: change the a, b, c to a', b', c' unit cell vectors
S_matrix = 1/2*[1 1 0; -1 1 0; 0 0 2];

% Step 3: change the coordinates according to new basis
coords = coords * S_matrix';

% Step 4: consider the symmetries and include quarter mirror images
coords_add = coords;
id = (coords(:,1) >= 0) & (coords(:,2) >= 0);
coords_add(id,1:2) = coords(id,1:2) - .5;

id = (coords(:,1) <= 0) & (coords(:,2) <= 0);
coords_add(id,1:2) = coords(id,1:2) + .5;

id = (coords(:,1) > 0) & (coords(:,2) < 0);
coords_add(id,1) = coords(id,1) - .5; coords_add(id,2) = coords(id,2) + .5;

id = (coords(:,1) < 0) & (coords(:,2) > 0);
coords_add(id,1) = coords(id,1) + .5; coords_add(id,2) = coords(id,2) - .5;

elems(Natom+1:2*Natom) = elems;

coords = [coords;coords_add];

% Step 5: translate the center from (1/2,1/2,0) to (0,0,0)
coords(:,1:2) = coords(:,1:2) + 0.5;

% Transforming back to initial fractional and then to cartesian.
trans_mat = S_matrix'^-1*trans_mat;
coords = coords * trans_mat;


end