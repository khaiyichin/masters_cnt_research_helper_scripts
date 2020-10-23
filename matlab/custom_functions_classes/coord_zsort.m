function out = coord_zsort(in)

[val,ind] = sort(in(:,3));
out = in(ind,:);
