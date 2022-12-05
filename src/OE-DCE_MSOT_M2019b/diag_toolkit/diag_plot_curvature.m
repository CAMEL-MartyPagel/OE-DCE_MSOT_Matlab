function diag_plot_curvature(fn)

load('-mat',fn);    

%%
figure;
subplot(2,2,1);
imagesc(max(BP,[],3));
axis image;
%     colormap gray;
set(gca,'XTick',xtick);
set(gca,'XTickLabel',num2str(xticklabel','%.1f'));
set(gca,'YTick',xtick);
set(gca,'YTickLabel',num2str(xticklabel','%.1f'));
xlabel('distance from transducer [mm]');
title(['MIP - backprojection (sensors ' num2str(sensors(1)) '-' num2str(sensors(end)) ')']);


%%
subplot(2,2,2);
imagesc(mip');
% mark peaks with lines
for i = 1:numel(imax)
    line([1 1]*imax(i),[0 10],'Color','white','LineStyle','-');
    line([1 1]*imax(i),[size(mip,2)-10 size(mip,2)],'Color','white','LineStyle','-');
end
set(gca,'YTick',ztick);
set(gca,'YTickLabel',num2str(zticklabel'));
ylabel('elevation [mm]');
set(gca,'XTick',xtick);
set(gca,'XTickLabel',num2str(xticklabel','%.1f'));
xlabel('distance from transducer [mm]');
%     grid on;
for j = 1:numel(ztick)
    line([0 size(mip,1)],[ztick(j) ztick(j)],'Color',[0.5 0.5 0.5],'LineStyle',':');
end
title('MIP along x');

%%
subplot(2,2,3);
imagesc(normsph');


set(gca,'YTick',ztick);
set(gca,'YTickLabel',num2str(zticklabel'));
ylabel('elevation [mm]');
xlabel([num2str(numpeaks) ' spheres']);
set(gca,'XTick',[]);
title('Normalised Elevational Field');
%     axis image;


%%
subplot(2,2,4);
hold off;
% add simulations
if ~isempty(simfile)
   load('-mat',simfile);
   area(simx*1e3+40,max([fwhm_min*zres; fwhm_max*zres]),'FaceColor','green','EdgeColor','none','DisplayName','Tolerance');
   hold on;
   area(simx*1e3+40,min([fwhm_min*zres; fwhm_max*zres]),'FaceColor','white','EdgeColor','none','DisplayName',' ');
   plot(simx*1e3+40,fwhm_opt*zres,'Color',[0 0.5 0],'DisplayName','Specification');
%        area(simx*1e3+40,fwhm_max*zres);

end

% Imasonic simualtion
try
    if (~isempty(imasim))
        IMsim = importdata(imasim);
        xres = 0.35;
        zres = 0.1;
        posx = 0:xres:xres*(size(IMsim,1)-1);
        posz = 0:zres:zres*(size(IMsim,2)-1);
         for d = 1:size(IMsim,1)
            vec = IMsim(d,:);
            vec = vec ./ max(vec);
            fwhm_i(d) = nnz(vec > 0.5);
         end
        plot(posx,fwhm_i*zres/2,'--','Color','red','DisplayName','Expectation');
        hold on;
    end
catch ex
        warning(['error loading Imasonic simulations: ' ex.message]);
end

plot(fwhm_dist*1e3,fwhm*zres,'.','Color','black','DisplayName','Measurements');
line([limits(1)*1e3 limits(2)*1e3],[1.2 1.2],'Color','black','LineStyle',':','DisplayName','Field Width 1.2mm');
ylim([0 max(fwhm*zres)])
xlim(limits*1e3);
set(gca,'XTick',xticklabel);
set(gca,'XTickLabel',num2str(xticklabel','%.1f'));
xlabel('distance from transducer [mm]');
ylabel('field width [mm]');
title('Field width');
legend('show','location','NorthWest');

