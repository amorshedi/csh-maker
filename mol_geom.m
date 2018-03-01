% gets coords and 2 or 3 indices and gives back bond length or angle
%the indices are vmd indieces

idx = [314 305 313] + 1;

if length (idx)==2
disp(['bond length is: ' num2str(norm(coords(idx(1),:)-coords(idx(2),:)))])
else

vec1 = coords(idx(1),:)-coords(idx(2),:);
vec2 = coords(idx(3),:)-coords(idx(2),:);

acosd(dot(vec1,vec2)/norm(vec1)/norm(vec2))

end