function Animat2_make_mat_file_allTrials
% Important variables to set: 
% weighted: should the average IIT values be weighted by the number of
% occurance of their respective state

clear all

Elem = 0:7;
%trialnum = [3 4 15 24 27 35 37 45 49 50 51 53 56 59 67 72 76 81 82 83 84 89 90 95 100 103 104 118 130 131 132 136 145 147 157 163 170 175 182 186 187 191 193 9 197];
%trialnum = [0 3 4 5 6 10 18 27 28 29 31 32 35 38 47 48 49 52 57 58 61 62 64 65 66 68 70 72 74 76 79 83 85 91 92 96 100 110 111 118 125 129 132 133 134 135 141 143 150 152 153 157 158 162 164 167 171 174 180 185 186 192 199];
%trialnum = [0 2 4 7 8 10 12 19 21 23 33 40 43 44 45 48 49 50 52 54 61 70 83 86 87 88 93 109 111 112 115 117 123 126 132 138 147 148 152 159 162 167 170 172 182 191 193 197];
%trialnum = [24 36 47 78 82 92 102 124 137 141 156 163 168 206 215 232 260 264 292 331 346 349 363 390]; %128 and 126
%trialnum = [24 36 47 92 102 137 141 156 163 215 232 264 292 331 346 349 363 390];
%trialnum = [7    10    19    23    45    48    50    54    55    93   109   111   117   126   138   152   162   170   172   182   191   193   197   206   213 220   240   258   285   330   344   368   372   399];
   
trialnum = 0:400;
numtrials = numel(trialnum);
totsteps = 60000-1;
step = 512;
range = 0:step:totsteps;
cond = 'c2a4_change_c23a14';

plotflag = 1;
weighted = 1;
DPath = '~/dev/animats/results/work_';
%DPath = '/Volumes/Macintosh HD 2/Simulations/Arend_XCodeAnimat/temporalSpatialIntegrationLite/work_';
path = strcat(DPath, cond, '/trial');

MaxFitness = 128;

ElemUsed_temp = zeros(1, numel(Elem));
FitPhiCorr = [];

count = 0;
evaluatedTrials = [];
for t = 1:numtrials
    Foldername = strcat(cond, '_trial', int2str(trialnum(t)));
    if exist(Foldername,'dir') == 7 
        count = count+1;
        evaluatedTrials = [evaluatedTrials t-1];
        for i = 1:length(range)
            %------------- get Fitness from Animat files---------------------------
            docname2 = strcat(path, int2str(trialnum(t)), '_', int2str(range(i)), '_KOdata.txt');

            Fitness = load(docname2);
            Fitness_level(count,i) = Fitness(1);  
    
            %------------- get rest from Phi calculation --------------------------
            % only go to folder here, because I need to access functions from the
            % main folder then, so I need to go out of the folder after loading
            cd(Foldername)
            Animat_gen = range(i);
            FilenameA = strcat(cond, '_animat_', int2str(trialnum(t)), '_', int2str(Animat_gen), '.mat');
            if exist(FilenameA,'file') == 2 
                load(FilenameA)
                cd ..
                
                num_lifeStates_temp = numel(PhiMip);
                % initialize 
                num_elements_complex = zeros(num_lifeStates_temp,1);
                if ~isempty(main_complex)
                    if num_lifeStates_temp == 1
                        num_elements_complex = numel(main_complex);
                    else
                        for j = 1:num_lifeStates_temp
                            % This is because the way python saves the
                            % main_complex, if all complexes have the same size
                            % it will be a row vector with several columns,
                            % while otherwise it's just a cell column vector
                            if size(main_complex, 1) > 1
                                num_elements_complex(j) = numel(main_complex(j,:));
                            else
                                num_elements_complex(j) = numel(main_complex{:,j});
                            end
                        end
                    end
                end
                
                if weighted == 1
                    % weighted by probability of occurance
                    p_LifeStates = double(cell2mat(LifeStates(:,2)));
                    p_LifeStates = p_LifeStates./sum(p_LifeStates);
                else
                    % just average
                    p_LifeStates = ones(length(PhiMip),1)./length(PhiMip);
                end
                
                all_used_nodes{count,i} = used_nodes;
                
                % double is necessary because output from Python is int64
                % and that can't be multiplied with doubles apparently
                big_phi_whole(count,i) = sum(p_LifeStates .* double(ssphi_whole_system)'); 
                big_phi_mip(count,i) = sum(p_LifeStates .* double(PhiMip)');  
                
                num_lifeStates(count,i) = length(PhiMip);
                num_connections(count,i) = number_of_connections;
                num_concepts_mip(count,i) = sum(p_LifeStates .* double(mip_num_concepts)');
                num_concepts_whole(count,i) = sum(p_LifeStates .* double(whole_num_concepts)');
                
                strongly_connected(count,i) = checkStrongCon(connectivity_matrix);
                
                size_main_complex(count,i) = sum(p_LifeStates .* num_elements_complex);
            else 
                cd ..
            end      
        end
    end
end    
%cd ..
%%
Fitness_Phi_Corr = corrcoef(mean(Fitness_level),mean(big_phi_mip))
MeanFitness = mean(Fitness_level,1);
MeanFitness = 100.*MeanFitness./MaxFitness;
%evaluatedTrials = trialnum;

File = strcat(cond, '_results');
save(File, 'weighted', 'evaluatedTrials','all_used_nodes','strongly_connected','big_phi_whole','big_phi_mip','Fitness_level','MeanFitness', 'num_lifeStates', 'num_connections', 'num_concepts_mip', 'num_concepts_whole', 'size_main_complex', 'cond', 'range');

if plotflag == 1
    figure
    subplot(3,1,1)
    hold on
    plot(range, MeanFitness)
    xlim([1, totsteps])
    
    subplot(3,1,2)
    hold on
    plot(range, mean(big_phi_whole,1), '-b')
    plot(range, mean(big_phi_mip,1), '-r')
    xlim([1, max(range)])

    subplot(3,1,3)
    hold on
    %plot(range, mean(num_lifeStates,1), '-b')
    plot(range, mean(num_connections,1), '-m')
    plot(range, mean(num_concepts_whole,1), '-b')
    plot(range, mean(num_concepts_mip,1), '-k')
    plot(range, mean(size_main_complex,1), '-g')
    xlim([1, max(range)])
end
end

function PotIrrComplex = checkStrongCon(connect_mat)
        J_sparse = sparse(connect_mat);
        [~,PotComplex] = graphconncomp(J_sparse);
        PotIrrComplex = length(unique(PotComplex))< length(connect_mat);
end
