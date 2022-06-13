function mcmcplot_density(chain,res,filepath_save,simulation_id)



names= res.names;

[nsimu,npar2]=size(chain);

inds = 1:npar2;

np  = length(inds);
ns1 = ceil(sqrt(np));
ns2 = round(sqrt(np));
for i=1:np
    h=subplot(ns1,ns2,i);
    hold on
    % [y,x]=density(chain(:,inds(i)),[],varargin{:});
    ksdensity_x = chain(:,inds(i));
    oldDir = pwd;
    cd('C:\Program Files\MATLAB\R2021b\toolbox\stats\stats\');
    h2 = str2func('ksdensity');
    [y,x] = ksdensity(ksdensity_x);
    cd(oldDir);

    y = y./max(y);

    plot(x,y,'-r')
    set(h,'ytick',[]);
    title(sprintf('%s',names{i}))
    
    box on
    grid on

    ylim([0,max(y)*1.2])
    yticks([0*max(y), .25*max(y), .5*max(y), .75*max(y), 1*max(y)])
    yticklabels({'0.00','0.25','0.50','0.75','1.00'})


    figure
    plot(x,y,'-r')
    set(h,'ytick',[]);
    title(sprintf('%s',names{i}))
    box on
    grid on
    ylim([0,max(y)*1.2])
    yticks([0*max(y), .25*max(y), .5*max(y), .75*max(y), 1*max(y)])
    yticklabels({'0.00','0.25','0.50','0.75','1.00'})
    saveas(gcf, [filepath_save  simulation_id  'density_' names{i}  '.fig'])
    saveas(gcf, [filepath_save simulation_id   'density_' names{i}  '.pdf'])
    saveas(gcf, [filepath_save simulation_id   'density_' names{i}  '.png'])
    print([filepath_save simulation_id 'density_' names{i}],'-depsc2','-painters')
    close

end


end