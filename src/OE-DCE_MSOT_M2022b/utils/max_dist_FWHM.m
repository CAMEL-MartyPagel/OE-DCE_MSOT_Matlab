function [D] = max_dist_FWHM(I,top_pix_num )

[x,y]=meshgrid(1:size(I,1),1:size(I,2));
Isorted=flip(sort(I(:)));
Imax=mean(Isorted(1:top_pix_num));
Ibin=(I>=0.5*Imax);   
x_selected=x(Ibin);
y_selected=y(Ibin); 

% figure;
% imagesc(Ibin)
% axis image

D=0;
for jj=1:length(x_selected)
    for ii=1:length(x_selected)
        D_test=sqrt((x_selected(jj)-x_selected(ii)).^2+(y_selected(jj)-y_selected(ii)).^2);               
        if D_test>=D;
            D=D_test;                    
        end                
    end
end


end

