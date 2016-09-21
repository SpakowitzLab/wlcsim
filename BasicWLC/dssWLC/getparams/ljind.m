function ind = ljind(l,j)
if (abs(j)<= l)
    ind = l^2+j+l+1;
else
    ind = 0;
end
end