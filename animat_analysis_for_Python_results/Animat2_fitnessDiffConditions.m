function Animat2_fitnessDiffConditions
clear all

global labels
labels = {};

options = [0, 1, 4];
taskA = {'c3a1'; 'c3a2';'c2a1'; 'c2a1'};
taskB = {'c14a23'; 'c36a45'; 'c23a14'};
%taskB = {'c23a14'; 'c36a45'; 'c14a23'};

for t = 1:length(taskA)
    for i = 1:length(taskB)
        for j = 1:length(options)
            option = options(j);
            condNum = i-1;
            if option == 0
                condition = char(strcat(taskB(i), '_36'));
                load(strcat(condition,'_results'));
            else
                condition = char(strcat(taskA(t), '_change_', taskB(i)));
                load(strcat(condition,'_results'));
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
            if t == 4
                ind2Sensors = zeros(1,length(evaluatedTrials));
                for l = 1:length(evaluatedTrials)
                    ind2Sensors(l) = length((intersect([0,1], all_used_nodes{l,59}))) == 2;
                end
            end
            if option == 1
                %indF = find(max(Fitness_level(:,55:59),[],2) == 128);
                %& (fitnessBeforeChange == 122)
                indF = find(max(big_phi_mip(:,55:59),[],2) > 0 & (fitnessBeforeChange >= 128) & (ind2Sensors' == 1));
                %Those with 2 nodes.
            %     indHU = [1     3     4     5     7     9    10    14    15    16    17    18    23    25    28    29    30    31    33    36    39    41    46    50    51    52    53    55    56    62];
            %     indF = intersect(indF, indHU);
                gcolor = [0 0 0];
            elseif option == 4
                %indF = find(max(Fitness_level(:,55:59),[],2) == 128);
                indF = find(max(big_phi_mip(:,55:59),[],2) == 0 & (fitnessBeforeChange >= 128) & (ind2Sensors' == 1));
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
            
            figure(200)
            subplot(length(taskA),length(taskB),1+condNum+(t-1)*length(taskB))
                function2use(inc, indX, rangeA, Fitness_level, gcolor, condition, [40,100])
                ylabel('Fitness (%)');

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