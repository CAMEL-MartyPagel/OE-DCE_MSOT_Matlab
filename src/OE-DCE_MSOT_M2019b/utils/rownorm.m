function y = rownorm(x);

y = x;
for j = 1:size(y,1),
    y(j,:) = y(j,:) - min(y(j,:));
    y(j,:) = y(j,:) ./ max(y(j,:));
end