
rng(50)
% clear %4 hours is needed i think...
% to_plot=0;
% 
% Fs=6;
% noiseThreshold=0.05;
% %load('/Volumes/Mac/RestoredMacBookPro/identification_paper/RAW/raw_data.mat');
% load('/Volumes/Mac/RestoredMacBookPro/identification_paper/HCTSA/QA_updated.mat');
% load('/Volumes/Mac/RestoredMacBookPro/identification_paper/RAW/raw_data_resp_sum.mat');
% Q=raw_data;
% SubjData=AllSubjData;
% SubjData=rmfield(SubjData,{'Data_sleep_HCTSA','Data_wake_HCTSA','mean_sleep','mean_wake','Data_LI_wake','Data_LI_sleep','Data_sleep','Data_wake','Pills'});

        int=5*60*Fs;
starts=randperm(7*21600,10000);

for perm=1:1000
start=starts(perm);
for i=1:size(Q,2)
    DataToUse= Q(i).wake;
    DataToUse(:,2)=- DataToUse(:,2);

    DataToUse_int=DataToUse(start:start+int-1,:);
[~,measureResults(i)]=NasalCycleParameters(DataToUse_int,Fs,noiseThreshold,'Mustache');
code_to_use=Q(i).Code;
tmp=strcmpi(code_to_use,{SubjData.Code});
demo(i)=SubjData(tmp);
end

 Cmat=table2array(struct2table(measureResults));
d=fields(demo);
m=fields(measureResults);
for j=2:size(d,1)
    TA=[demo.(d{j})];
    for ii=1:size(Cmat,2)
        u=unique(TA);
        u(isnan(u))=[];
        if size(u,2)>3

            x = Cmat(:,ii);
y = TA(:);

% base validity
keep = isfinite(x) & isfinite(y);

% remove outliers ±2.5 SD (on x and y)
k = 2.5;

mx = mean(x(keep),'omitnan'); sx = std(x(keep),'omitnan');
my = mean(y(keep),'omitnan'); sy = std(y(keep),'omitnan');

if isfinite(sx) && sx>0
    keep = keep & (abs(x - mx) <= k*sx);
end
if isfinite(sy) && sy>0
    keep = keep & (abs(y - my) <= k*sy);
end

% correlation on cleaned data
[R,p] = corr(x(keep), y(keep), 'rows','complete');
 
RR(perm,j)=R;
p_corr_perm(perm,j)=p;
fprintf(['significant correlaton between ' d{j} ' and ' m{ii} ': R = ' num2str(R) ' p = ' num2str(p) ' \n'])
            % plot  
             if p<0.001 && to_plot

                figure('Color','w');
                scatter(Cmat(:,ii), TA,60,'filled');
                hold on
                h = lsline;
set(h, 'LineWidth', 2.5, 'Color','k');  

                xlabel(m{ii}, 'Interpreter','none');
                ylabel(d{j}, 'Interpreter','none');
                title(sprintf('%s vs %s | R=%.3f, p=%.4g,', ...
                    d{j}, m{ii}, R, p), 'Interpreter','none');
            
            end
        elseif size(u,2)==2
             x = Cmat(:,ii);
    g = TA(:);  % numeric 0/1 expected

    % keep only valid rows + only groups 0/1
    keepBase = isfinite(x) & isfinite(g) & ismember(g,[0 1]);

    % --- remove outliers: ±2.5 SD WITHIN EACH GROUP ---
    k = 2.5;
    keep = false(size(x));

    for gg = [0 1]
        idx = keepBase & (g==gg);
        xv  = x(idx);

        mu = mean(xv,'omitnan');
        sd = std(xv,'omitnan');

        if ~isfinite(sd) || sd==0
            keep(idx) = true;  % nothing to remove
        else
            keep(idx) = abs(xv - mu) <= k*sd;
        end
    end

    x_clean = x(keep);
    g_clean_num = g(keep);

    % ttest2 on cleaned data
    [~,pp] = ttest2(x_clean(g_clean_num==0), x_clean(g_clean_num==1));
ppp(perm,j)=pp;
        fprintf('significant ttest2 between %s and %s (±%.1f SD removed): p = %.4g (n=%d)\n', ...
            d{j}, m{ii}, k, pp, numel(x_clean));
    if pp < 0.1 && to_plot

        % --- violin + dots ---
        grp_clean = categorical(g_clean_num, [0 1], {'0','1'});
        figure('Color','w'); hold on

        vals0 = x_clean(grp_clean=='0');
        vals1 = x_clean(grp_clean=='1');

 violin({vals0(:), vals1(:)}, {'0','1'});

% get violin patches
vp = findobj(gca, 'Type','Patch');

% example colors (one per group)
cols = [ 0.9 0.5 0.4; ...
    0.4 0.6 0.9];  % red-ish

for kk = 1:numel(vp)
    vp(kk).FaceColor = cols(kk,:);
    vp(kk).FaceAlpha = 0.6;
    vp(kk).EdgeColor = 'k';
end
        xpos = double(grp_clean);
        jit  = (rand(size(xpos))-0.5)*0.12;
        plot(xpos + jit, x_clean, 'k.', 'MarkerSize', 10);

        % medians
        % plot(1, median(x_clean(g_clean_num==0),'omitnan'), 'kd', 'MarkerFaceColor','k');

        % plot(2, median(x_clean(g_clean_num==1),'omitnan'), 'kd', 'MarkerFaceColor','k');
        xlabel(d{j}, 'Interpreter','none');
        ylabel(m{ii}, 'Interpreter','none');
        title(sprintf('%s by %s | ttest2 p=%.4g | ±%.2f SD | n0=%d n1=%d', ...
            m{ii}, d{j}, pp, k, sum(g_clean_num==0), sum(g_clean_num==1)), ...
            'Interpreter','none');   
        legend off

    end
        else
            disp(unique(TA))
        end

    end
end

end
