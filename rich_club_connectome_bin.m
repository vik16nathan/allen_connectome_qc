% load('demo_networks.mat')
% from a connectivity matrix we want to obtain:
% rich club: the rich club regime, which nodes are part of the rich club,
% frequency of each node being part of the rich club.
clear
%%add path for Brain Connectivity Toolbox
addpath("2019_03_03_BCT/");
output_dir='/data/chamal/projects/natvik/knox_qc_full_06232025/derivatives/rich_club/';
%%
knox_rgn_conn_matrices=["../derivatives/regionalized_connectomes/knox_conn_strength_old_bilateral.csv",...
    "../derivatives/regionalized_connectomes/knox_conn_strength_new_bilateral_1018.csv"];
knox_rgn_conn_suffixes=["old", "new_1018"];

for percentile_threshold=[70, 80, 90]
    for i=1:2 %%%run for old/new connectomes
        %%set same random seed
        rng(123); % Sets the seed to 123
        knox_rgn_conn_matrix_original=readmatrix(knox_rgn_conn_matrices(i));
        suffix=char(knox_rgn_conn_suffixes(i));
        %get rid of indices
        sc=knox_rgn_conn_matrix_original;
        sc=sc(2:size(sc,1),2:size(sc,2));
        n = length(sc);
        sc = sc.*~eye(n);
        numeric_thresh=prctile(sc, percentile_threshold, 'all');
        sc(sc > numeric_thresh) = 1; % binarize using 80% threshold (keep top 20% conns)
        sc(sc <= numeric_thresh) = 0; % binarize using 80% threshold (keep top 20% conns)
        
        %%
        d = distance_bin(sc);
        [ind,outd,degree] = degrees_dir(sc);
        
        %%
        %--------------------------------------------------------------------------
        % find rich club in a structural network
        %--------------------------------------------------------------------------
        nswap = 2^6;    % no. of flips/edge
        nrand = 1000;     % no. of randomized networks
        
        k = sum(sc,2);                  % degree of each node
        %kmax = 118;                      % max degree
        [r,nk,ek] = rich_club_bd(sc,max(degree));   %max(degree) = 190
        
        %%
        
        rrand = zeros(nrand,max(degree));
        nkrand = zeros(nrand,max(degree));
        ekrand = zeros(nrand,max(degree));
        
        for irand = 1:nrand
            b = randmio_dir(sc,nswap);
            [rrand(irand,:),nkrand(irand,:),ekrand(irand,:)] = rich_club_bd(b,max(degree));
            
            fprintf('rand %i of %i done\n',irand,nrand);
        end
        %% replace inf or nan with 0
        rrand(isinf(rrand)|isnan(rrand)) = 0;
        r(isinf(r)|isnan(r)) = 0;
        
        %%
        kmax=max(degree);
        pval = zeros(1,kmax);
        pval_norm = zeros(1,kmax);
        
        r_norm = zeros(1,kmax);
        for ik = 1:kmax
            r_norm(ik) = r(ik)/mean(rrand(:,ik));
        end
        for ik = 1:kmax
            pval(ik) = (sum(rrand(:,ik)>=r(ik))/nrand) ;
        end
        r_norm(isinf(r_norm)|isnan(r_norm)) = 0;
        
        %% --------------------------------------------------------------------------
        % Final plot for rich club
        %--------------------------------------------------------------------------
        r_norm_sig=NaN(size(r_norm));
        r_norm_sig(find(pval<=0.05))=r_norm(find(pval<=0.05)); 
        figure;
        
        plot(r_norm_sig','MarkerFaceColor',[1 1 1],'LineStyle','-','LineWidth',2,'Marker','.','MarkerSize',35,'Color','#e37878');
        hold on
        plot(r','MarkerFaceColor',[1 1 1],'LineStyle','-','LineWidth',2,'Marker','.','MarkerSize',15,'Color','#474747');
        hold on
        plot(mean(rrand)','MarkerFaceColor',[1 1 1],'LineStyle','-','LineWidth',2,'Marker','.','MarkerSize',15,'Color','#8c8b8b');
        
        plot(r_norm,'MarkerFaceColor',[1 1 1],'LineStyle','-','LineWidth',2,'Marker','.','MarkerSize',10,'Color','#cf0606');
        yline(1);
        axis([0 kmax min(r_norm)-.1 max(r_norm)+.2])
        legend('pval<0.05','rich club','rich club rand','normalized rich club')
        box off
        title('Rich Club')
        ax=gca;
        ax.LineWidth=1.5;
        set(gcf,'color','w');
        xlabel('Degree')
        saveas(gcf, [output_dir,'rich_club_',suffix,'_pct',char(percentile_threshold),'.png']);
        hold off
        %% -------------------------------------------------------------------------
        % topological rich club regime: 
        %--------------------------------------------------------------------------
        degrez=1:length(r_norm_sig);
        sig_degrees=~isnan(r_norm_sig);
        sig_degrees2=degrez.*sig_degrees;
        sig_degrees2(sig_degrees2==0)=NaN;
        topo=mean(sig_degrees2,'omitmissing')+std(sig_degrees2,'omitmissing');
        fprintf('Topo %6.2f \n',topo);
        %% -------------------------------------------------------------------------
        % define frequency of which nodes are rich club
        %--------------------------------------------------------------------------
        degrees_loop=sig_degrees2(~isnan(sig_degrees2));
        richOrNot_mult=zeros(size(degrees_loop,2),size(sc,1));
        
        for j=1:size(degrees_loop,2)
            for i=1:size(sc,1)
                if degree(i)>=degrees_loop(j)
                    richOrNot_mult(j,i)=1;
                end
            end
        end
        richOrNot_mult_final=sum(richOrNot_mult,1);
        richOrNot_mult_final_ids_pctg=richOrNot_mult_final/max(richOrNot_mult_final(1,:));
        save([output_dir,'knox_conn_',suffix,'_pct',char(percentile_threshold), '_richOrNot_pctg.mat'] , 'richOrNot_mult_final_ids_pctg');
        
        
        %% -------------------------------------------------------------------------
        % define which nodes are rich club
        %--------------------------------------------------------------------------
        nodes_ids=1:1:size(sc,1);
        richOrNot=zeros(1,size(sc,1));
        for i=1:size(sc,1)
            if degree(i)>=topo
                richOrNot(i)=1;
            end
        end
        save([output_dir,'knox_conn_',suffix,'_pct',char(percentile_threshold), '_richOrNot_topo.mat'],'richOrNot');
    end
end
