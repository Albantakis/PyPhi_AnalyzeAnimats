function Animat2_plotData
clear all

global labels
labels = {};

Hamming = 0;
InfoMeasures = 0;
PlotHistagrams = 1;
xbins = 63:1:128

options = [0, 1, 4];
taskA = {'c2a4'};% 'c3a1'; 'c3a2'; } %;'c3a1'};
taskB = {'c14a23','c23a14'};% 'c36a45'; 
%taskB = {'c23a14'; 'c36a45'; 'c14a23'};

for i = 1:length(taskB)
    for j = 1:length(options)
        option = options(j);
        condNum = i-1;
        if option == 0
            condition = char(strcat(taskB(i), '_36'));
            load(strcat(condition,'_results'));
        else
            if length(taskA) > 1
                for k = 1:length(taskA)
                    conditiontemp{k} = char(strcat(taskA(k), '_change_', taskB(i)));                   
                    data{k} = load(strcat(conditiontemp{k},'_results'));                   
                end
                % specify dimension across which the .mat file Variables
                % should be concatenated
                dimIndex = ones(numel(fieldnames(data{1})),1);
                dimIndex([2,4]) = 2;
                %dimIndex([2,3]) = 2;
                newMatfile = concatenate_matfiles_from_struct(data, dimIndex);
                condition = newMatfile.cond;
                save(condition,'-struct','newMatfile')
                load(strcat(condition));
            else
                condition = char(strcat(taskA, '_change_', taskB(i)));
                load(strcat(condition,'_results'));
            end
        end
        %condition = 'c36a45_36';

        %function2use = @plotIndividualTraces;
        function2use = @plotVariable;
        %evaluatedTrials = evaluatedTrials(evaluatedTrials < 300);

        %% Data
        for l = 1:length(evaluatedTrials)
            labels = [labels; int2str(evaluatedTrials(l))];
        end
        fitnessBeforeChange = max(Fitness_level(:,55:59),[],2);
        ind2Sensors = ones(1,length(evaluatedTrials));

%         for l = 1:length(evaluatedTrials)
%             ind2Sensors(l) = length((intersect([0,1], all_used_nodes{l,59}))) == 2;
%         end
         
        if option == 1
            %indF = find(max(Fitness_level(:,55:59),[],2) == 128);
            %& (fitnessBeforeChange == 122)
            indF = find(max(big_phi_mip(:,55:59),[],2) > 0 & (fitnessBeforeChange >= 128) & (ind2Sensors' == 1));
            %indF = find(max(strongly_connected(:,55:59),[],2) > 0 & (fitnessBeforeChange >= 128) & (ind2Sensors' == 1));
            %Those with 2 nodes.
        %     indHU = [1     3     4     5     7     9    10    14    15    16    17    18    23    25    28    29    30    31    33    36    39    41    46    50    51    52    53    55    56    62];
        %     indF = intersect(indF, indHU);
            gcolor = [0 0 0];
        elseif option == 4
            %indF = find(max(Fitness_level(:,55:59),[],2) == 128);
            indF = find(max(big_phi_mip(:,55:59),[],2) == 0 & (fitnessBeforeChange >= 128) & (ind2Sensors' == 1));
            %indF = find(max(strongly_connected(:,55:59),[],2) == 0 & (fitnessBeforeChange >= 128) & (ind2Sensors' == 1));
            %Those with more than 2 nodes.
        %     indHU = [2     6     8    11    12    13    19    20    21    22    24    26    27    32    34    35    37    38    40    42    43    44    45    47    48    49    54    57    58    59    60    61];
        %     indF = intersect(indF, indHU);
            gcolor = [0 0 1];
        else
            indF = [1:size(Fitness_level, 1)];
            gcolor = [0.75 0.75 0.5];
            %gcolor = [0 0.5 1];
        end

        totsteps = 60000-1;

        inc = reshape(indF, 1, []); %[1:50];
        inc = inc(inc <= length(evaluatedTrials))
        length(inc)

        MaxFitness = 128;
        % MeanFitness = mean(Fitness_level(inc,:),1);
        % MeanFitness = 100.*MeanFitness./MaxFitness;
        Fitness_level = 100.*Fitness_level./MaxFitness;

        %% Fitness and Phi
        rangeB = eval('range');
        rangeB = rangeB(1,:);
        indX = 1:length(rangeB);
        rangeA = rangeB(1,indX)./10000;
        figure(34)

        subplot(6,4,1+condNum)
            function2use(inc, indX, rangeA, Fitness_level, gcolor, condition, [40,100])
            ylabel('Fitness (%)');
        subplot(6,4,5+condNum)
            function2use(inc, indX, rangeA, num_concepts_whole, gcolor, condition, [0,5])
            ylabel('#concepts');
        subplot(6,4,9+condNum)
            function2use(inc, indX, rangeB, big_phi_whole, gcolor, condition, [0,1])
            ylabel('$\sum \varphi^{\rm Max}$','Interpreter','latex','FontSize',12);

        subplot(6,4,13+condNum)
            function2use(inc, indX, rangeA, size_main_complex, gcolor, condition, [0,1.5])
            ylabel('#elements');
        subplot(6,4,17+condNum)
            function2use(inc, indX, rangeA, num_concepts_mip, gcolor, condition, [0,2])
            ylabel('#concepts');
        subplot(6,4,21+condNum)
            function2use(inc, indX, rangeB, big_phi_mip, gcolor, condition, [0,0.4])
            ylabel('$\Phi^{\rm Max}$','Interpreter','latex','FontSize',12);
            xlabel('#Generations');

        %% I_SMMI and I_Pred
        if InfoMeasures == 1
            cd ..
            cd('I_sen_pred')
            load(strcat(condition,'_ISenMot'));
            I_SMMI = results.ISenMot;
            HMot = results.HMot;
            HSen = results.HSen;

            load(strcat(condition,'_IPred'));
            I_Pred = results.IPred;
            H_States = results.HStates;

            figure(102)
            subplot(6,4,1+condNum)
                function2use(inc, indX, rangeA, Fitness_level, gcolor, condition, [40,100])
                ylabel('Fitness (%)');
            subplot(6,4,5+condNum)
                function2use(inc, indX, rangeA, I_SMMI, gcolor, condition, [0,1])
                ylabel('I_{SMMI}')
            subplot(6,4,9+condNum)
                function2use(inc, indX, rangeA, HSen, gcolor, condition, [0,2])
                ylabel('H_{Sen}')   
            subplot(6,4,13+condNum)
                function2use(inc, indX, rangeA, HMot, gcolor, condition, [0,2])
                ylabel('H_{Mot}')
            subplot(6,4,17+condNum)
                function2use(inc, indX, rangeA, I_Pred, gcolor, condition, [0,4])
                ylabel('I_{Pred}')
            subplot(6,4,21+condNum)
                function2use(inc, indX, rangeB, H_States, gcolor, condition, [0,4])
                ylabel('H_{State}')
                xlabel('#Generations');
            cd ..
            cd('weighted')
        end
        
        if Hamming == 1 && option ~= 0
            cd ..
            cd(strcat('Hamming/', taskA{1}))
%             load(strcat('HammingDistCM_', condition, '_BC'));           
%             UniqueCMs = unique(WhichUniqueCM(inc));
%             numUniqueCMs = numel(UniqueCMs)
%             for k = UniqueCMs
%                 %WhichUniqueCM gives unique index#s, the index# plugged
%                 %into UniquecmsInd gives me the generation number that fits
%                 %the example cm, +1 s are because the indices are python
%                 view(biograph(reshape(Allncm2D(1+UniquecmsInd(k+1),:),8,[])'))
%             end
            HCM = load(strcat('HammingDistCM_', condition, '_118_prunedFalse'));
            HTPM = load(strcat('HammingDistTPM_', condition, '_118_prunedFalse'));
            
            HCM_prunedTrue = load(strcat('HammingDistCM_', condition, '_118_prunedTrue'));
            HTPM_prunedTrue = load(strcat('HammingDistTPM_', condition, '_118_prunedTrue'));
            
            figure(103)
            subplot(5,4,1+condNum)
                function2use(inc, indX, rangeA, Fitness_level, gcolor, condition, [40,100])
                ylabel('Fitness (%)');
            subplot(5,4,5+condNum)
                function2use(inc, indX, rangeA, double(HCM.HammingDist), gcolor, condition, [0,10])
                ylabel('HDist CM')
            subplot(5,4,9+condNum)   
                function2use(inc, indX, rangeA, double(HCM_prunedTrue.HammingDist), gcolor, condition, [0,8])
                ylabel('HDist pruned CM')
            subplot(5,4,13+condNum)
                function2use(inc, indX, rangeA, double(HTPM.HammingDist), gcolor, condition, [0,500])
                ylabel('Hist TPM')
            subplot(5,4,17+condNum)
                function2use(inc, indX, rangeA, double(HTPM_prunedTrue.HammingDist), gcolor, condition, [0,300])
                ylabel('HDist pruned TPM')
            
            cd ..
            cd ..
            cd('weighted')
        end
        
        
        if PlotHistagrams == 1 && option ~= 0
            rangeToPlot = [60, 70, 80, 90, 100, 118];
            figure(104)
            for l = 1:6
            subplot(6,4,(l-1)*4+condNum+1)
                hold on
                hist(Fitness_level(inc,rangeToPlot(l)));
                h = findobj(gca, 'Type', 'patch');
                h = h(1);
                if option == 1
                    h.FaceColor = [gcolor];
                elseif option == 4
                    h.EdgeColor = [gcolor];
                    h.FaceColor = 'none';
                end
                YLIM = [0 15];
                set(gca,'ylim',YLIM);
                XLIM = [min(xbins)-0.5 max(xbins)+0.5];
                set(gca,'xlim',XLIM);
            end
        end
    end
end
end

function plotVariable(inc, indX, rangeX, variable, gcolor, condition, yLimits)
  hold on
  if numel(inc) == 1
    SEM = zeros(1, length(variable));
  else  
    SEM = std(variable(inc,:),1)/sqrt(length(inc));
  end  
    shadedErrorBar(rangeX(indX), mean(variable(inc,indX),1), SEM(indX), {'color', gcolor,'DisplayName',condition,'LineWidth', 1},0);
    xlim([0, max(rangeX)])
    ylim(yLimits);
end   

function plotIndividualTraces(inc, indX, rangeX, variable, gcolor, condition, yLimits)
  global labels
  hold on
  p = plot(rangeX(indX), variable(inc,indX), 'LineWidth', 1);
  set(p, {'DisplayName'}, labels(inc));
  xlim([0, max(rangeX)])
  ylim(yLimits);
end  