function cax = autothres(img)

    [ne cent] = hist(img(:),100);
    cum = ne*triu(ones(100,100),0); cum = cum ./ max(cum(:));
    crophigh = find(cum > 0.999); % find elements exceeding % threshold
%     croplow = find(cum < 0.01); if isempty(croplow), croplow = 1; end
%     cax = [cent(croplow(end)) cent(crophigh(1))]; % use center as thres
    cax = [0 cent(crophigh(1))];

