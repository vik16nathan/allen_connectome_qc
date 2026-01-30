%% -------------------------------------------------------------------------
% find communities using Louvain algorithm
%--------------------------------------------------------------------------
%%load whole-brain REGIONALIZED structural connectivity .csv 
cd(fileparts(mfilename('fullpath')))
addpath("2019_03_03_BCT/");
output_dir='../derivatives/community_louvain/';
rich_club_dir='../derivatives/rich_club/';
%%
knox_rgn_conn_matrices={'../derivatives/regionalized_connectomes/knox_conn_strength_old_bilateral.csv',...
    '../derivatives/regionalized_connectomes/knox_conn_strength_new_bilateral.csv',...
    '../derivatives/regionalized_connectomes/knox_conn_density_old_bilateral.csv',...
    '../derivatives/regionalized_connectomes/knox_conn_density_new_bilateral.csv'};

knox_rgn_conn_suffixes={'old_strength', 'new_strength', 'old_density', 'new_density'};
%Repeat analysis for old vs. new connectomes
%For consistency, assign gamma with the most stable clustering (largest cluster of high mutual information)
%as the same (gamma index=21) b/t old and new connectomes
%%right now, assuming it's the same @ 21 for both connectomes
%stick this in the supplement


%Binary threshold for minimal agreement across gammas to visualize
agreement_thresh=0.5; 
%%
%%for percentile_threshold=[70, 80, 90]
for percentile_threshold=[80]
    for i=1:numel(knox_rgn_conn_matrices)
        %%binarize structural connectivity
        knox_rgn_conn_matrix_original=readmatrix(knox_rgn_conn_matrices{i});
        knox_rgn_conn_names=knox_rgn_conn_matrix_original(1,2:size(knox_rgn_conn_matrix_original, 2));
        suffix=char(knox_rgn_conn_suffixes{i});
        %get rid of indices
        sc=knox_rgn_conn_matrix_original;
        sc=sc(2:size(sc,1),2:size(sc,2));
        n = length(sc);
        sc = sc.*~eye(n);
        numeric_thresh=prctile(sc, percentile_threshold, 'all');
        sc(sc > numeric_thresh) = 1; % binarize using % threshold (keep top 1 - % conns)
        sc(sc <= numeric_thresh) = 0; % binarize using % threshold (keep top 1 - % conns)
        %%
        nodes_ids=1:1:size(sc,1);
        gamvals = 0.3:0.05:3;
        ngam = length(gamvals);
        
        nreps = 100;
        ci = zeros(n,ngam,nreps);   % community assignments
        q = zeros(ngam,nreps);      % modularity value
        ciu = zeros(n,ngam);        % consensus communities
        
        % parpool;
        for igam=1:ngam
            gam = gamvals(igam);
            B = sc;
            parfor irep = 1:nreps
                    [ci(:,igam,irep),q(igam,irep)] = community_louvain(B,gamvals(igam),[],[]);
            end
            citemp = squeeze(ci(:,igam,:)); 
            if igam==15 %%make plot at a single igam to represent community assignments
               figure
               imagesc(citemp)
               colormap(jet)
               xlabel('Repetitions')
               ylabel('Nodes')
               title('Community assignments (one gamma value)')
               colorbar
               set(gcf,'color','w');
            end
        %%% community assignments for a single gamma value (nodes x repetitions)
            ag = agreement(citemp) / nreps; % agreement (probability)
            
            cinull = zeros(size(citemp));
        % get null agreement
            parfor irep = 1:nreps
                cinull(:,irep) = citemp(randperm(n),irep);
            end
            
        
            agnull = agreement(cinull) / nreps;
            tau = mean(agnull(:));
            %tau = mean(ag(:));
        
            % get consensus (see Lancichinetti & Fortunato (2012))
            ciu(:,igam) = consensus_und(ag,tau,10);
        end
        % plot consensus community assignments
        figure;
        imagesc(ciu); 
        xlabel('gamma values')
        ylabel('nodes')
        colormap(jet)
        title('Community Assignments (one for each gamma)')
        colorbar;
        set(gcf,'color','w');
        %% -------------------------------------------------------------------------
        % Consensus partition using adjusted mutual information
        %--------------------------------------------------------------------------
        rand_input_mat_path=[output_dir,'rand_input_',suffix,'.mat'];
        save(rand_input_mat_path,'ciu');
        %run this in a separate shell
        system(['python3 call_randind.py --matfile ', rand_input_mat_path, ' --flag 1 --suffix ', suffix]);
        %%
        x=load(['../derivatives/community_louvain/rand_output_',suffix,'.mat']);
        %%
        ami_thr=0.9:0.005:1;
        ami_matrix=x.a;
        %find consensus clusters within gamma x gamma mutual information matrix,
        %at a variety of different mutual information thresholds
        ciu_ami = zeros(size(ciu,2),length(ami_thr)); 
        
        for i=1:length(ami_thr)
            % agfinal(agfinal<=0.80) = ami_thr(1);
            % consensus_und(ag,tau,10);
            ciu_ami(:,i) = consensus_und(ami_matrix,ami_thr(i),10);
        end
    
        %%
        x.a(x.a<=0.9) = 0;
        %% -------------------------------------------------------------------------
        % Plot raw adjusted mutual information
        %--------------------------------------------------------------------------
        figure;
        imagesc(x.a)
        colormap(hot)
        xlabel('Partitions')
        ylabel('Partitions')
        title('Adjusted Mutual Information')
        colorbar;
        set(gcf,'color','w');
        axis square;
        caxis([0, 1]);
        
        
        %% -------------------------------------------------------------------------
        % Threshold adjusted mutual information
        %--------------------------------------------------------------------------
        agfinal=agreement(ciu_ami) /length(ami_thr);
        agfinal=agfinal+eye(length(gamvals));
        agfinal(agfinal<=agreement_thresh) = 0;
        % agfinal = agfinal.*eye(length(agfinal));
        figure
        imagesc(agfinal)
        colormap(flip(gray))
        xlabel('Gamma Index for Partition')
        ylabel('Gamma Index')
        title('Agreement Between Consensus Partitions')
        colorbar;
        set(gcf,'color','w');
        axis square;
        saveas(gcf,[output_dir,'adjusted_mutual_inf_',suffix,'_pct',num2str(percentile_threshold),'.png']);
        
        %% -------------------------------------------------------------------------
        % Select gamma and particion data
        %--------------------------------------------------------------------------
        % select interval of gamma
        % caculate dice coeff betweeen partitoon in gamma interval
        % select single gamma
        
        % choose range between gamma=0.65,0.7,0.75,0.8,0.85,0.9
        % for a single gamma value:
    
        
        B = sc;
        nreps=100;
        n=length(B);
        com_asigns=zeros(length(B),nreps);
        q=zeros(length(B),nreps);
        
        gamma_chosen=21; %%%%index of "center" of largest cluster of mutual information scores (see above)
        parfor irep = 1:nreps
            [com_asigns(:,irep),q(irep)] = community_louvain(B,gamvals(gamma_chosen),[],[]);
        end
        
        ag = agreement(com_asigns) / nreps; % agreement (probability)
            
        cinull = zeros(size(com_asigns));
        % get null agreement
        parfor irep = 1:nreps
            cinull(:,irep) = com_asigns(randperm(n),irep);
        end
        agnull = agreement(cinull) / nreps;
        tau = mean(agnull(:));
        %tau = mean(ag(:));
        
        % get consensus (see Lancichinetti & Fortunato (2012))
        final_clust(:) = consensus_und(ag,tau,10); 
        % save final cluster assignments
        metadata.('comm_assign')=final_clust';
        %% -------------------------------------------------------------------------
        % Create csv for gephi (source-target columns)
        %--------------------------------------------------------------------------
        % save for gephi 
        [trow,tcol]=find(sc); %%anything  > 0
        nodes=[tcol,trow];
        Source=tcol;
        Target=trow;
        nodes_table=table(Source,Target);
        metadata.names=knox_rgn_conn_names;
        %metadata.Id=unique(Source);
        metadata.Id=1:n;
        writetable(nodes_table,[output_dir, 'gephi_source_target_',suffix,'_pct',num2str(percentile_threshold),'.csv']);
        
        %% define participantion coeff
        P_out=participation_coef(sc,metadata.comm_assign,1);
        P_in=participation_coef(sc,metadata.comm_assign,2);
        metadata.participation_out=P_out;
        metadata.participation_in=P_in;
        v=P_in;
        max_P=1-(1/(max(metadata.comm_assign)));
        metadata.part_in_ptg=(P_in.*1)/max_P;
        metadata.part_out_ptg=(P_out.*1)/max_P;
        
        %% -------------------------------------------------------------------------
        % Save updated csv with all new metadata + rich club
        %--------------------------------------------------------------------------
        rich_club_pct=load([rich_club_dir,'knox_conn_',suffix, '_pct',num2str(percentile_threshold), '_richOrNot_pctg.mat']);
        metadata.rich_club_pct=rich_club_pct.richOrNot_mult_final_ids_pctg;
    
        rich_club_topo=load([rich_club_dir,'knox_conn_',suffix,'_pct',num2str(percentile_threshold), '_richOrNot_topo.mat']);
        metadata.rich_club_topo=rich_club_topo.richOrNot;
        % include rich club assignments from rich_club_vikram.m
        names = fieldnames(metadata);
        T = table();
        
        for k = 1:numel(names)
            T.(names{k}) = metadata.(names{k})(:);  % ensure column vector
        end
        
        % Write to CSV
        writetable(T, [output_dir,'metadata_community_rich_club_',suffix,'_pct',num2str(percentile_threshold),'.csv']);  % change filename as needed
    end
end
