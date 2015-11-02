function Animat2_plotData_individual
clear all

global labels
labels = {};

DPath = '/Users/larissa_alb/dev/animats/results/work_';
%DPath = '/Volumes/Macintosh HD 2/Simulations/Arend_XCodeAnimat/temporalSpatialIntegrationLite/work_';

plot_network = 1

option = 1;
condition = 'c3a1_change_c36a45';
%condition = 'c36a45_36';

function2use = @plotIndividualTraces;
I_sen_mot = 0;
%% Data
load(strcat(condition,'_results'));

for l = 1:length(evaluatedTrials)
    labels = [labels; int2str(evaluatedTrials(l))];
end
fitnessBeforeChange = max(Fitness_level(:,59),[],2);

if option == 1
    %indF = find(max(Fitness_level(:,55:59),[],2) == 128);
    %& (fitnessBeforeChange == 122)
    indF = find(max(big_phi_mip(:,55:59),[],2) > 0 & (fitnessBeforeChange >= 126));% & max(BigPhiMip(:,60:64),[],2) == 0) % find(mean(Fitness_level(:,end-10:end),2) > FA(i) & mean(Fitness_level(:,end-10:end),2) <= FA(i+1))
    %Those with 2 nodes.
%     indHU = [1     3     4     5     7     9    10    14    15    16    17    18    23    25    28    29    30    31    33    36    39    41    46    50    51    52    53    55    56    62];
%     indF = intersect(indF, indHU);
    gcolor = [0 0 0];
elseif option == 4
    %indF = find(max(Fitness_level(:,55:59),[],2) == 128);
    indF = find(max(big_phi_mip(:,55:59),[],2) == 0 & (fitnessBeforeChange >= 126));% & max(BigPhiMip(:,60:64),[],2) == 0) % find(mean(Fitness_level(:,end-10:end),2) > FA(i) & mean(Fitness_level(:,end-10:end),2) <= FA(i+1))
    %Those with more than 2 nodes.
%     indHU = [2     6     8    11    12    13    19    20    21    22    24    26    27    32    34    35    37    38    40    42    43    44    45    47    48    49    54    57    58    59    60    61];
%     indF = intersect(indF, indHU);
    gcolor = [0 0 1];
else
    indF = [1:size(Fitness_level, 1)];
    %gcolor = [0.75 0.75 0.5];
    gcolor = [0 0.5 1];
end

totsteps = 60000-1;

length(indF)

inc = [17, 19, 24, 35, 43];%indF; %[1:50];
inc = reshape(inc, 1, []) %[1:50];

MaxFitness = 128;
% MeanFitness = mean(Fitness_level(inc,:),1);
% MeanFitness = 100.*MeanFitness./MaxFitness;
%Fitness_level = 100.*Fitness_level./MaxFitness;

%% Fitness and Phi
rangeB = eval('range');
indX = 1:length(rangeB);
rangeA = rangeB(indX)./10000;
nrows = 3;
ncols = 2;
for i = inc
    figure(1000+i)
    subplot(nrows,ncols,1)
        plotVariable(inc, indX, rangeA, Fitness_level, [.3 .3 .3] , condition, [60,128])
        function2use(i, indX, rangeA, Fitness_level, gcolor, condition, [60,128])
        ylabel('Fitness (%)');
    subplot(nrows,ncols,3)
        plotVariable(inc, indX, rangeA, num_concepts_whole, [.3 .3 .3] , condition, [0,7])
        function2use(i, indX, rangeA, num_concepts_whole, gcolor, condition, [0,7])
        ylabel('#concepts');
    subplot(nrows,ncols,5)
        plotVariable(inc, indX, rangeB, big_phi_whole, [.3 .3 .3] , condition, [0,1])
        function2use(i, indX, rangeB, big_phi_whole, gcolor, condition, [0,1])
        ylabel('$\sum \varphi^{\rm Max}$','Interpreter','latex','FontSize',12);
    subplot(nrows,ncols,2)
        plotVariable(inc, indX, rangeA, size_main_complex, [.3 .3 .3] , condition, [0,2.5])
        function2use(i, indX, rangeA, size_main_complex, gcolor, condition, [0,2.5])
        ylabel('#elements');
    subplot(nrows,ncols,4)
        plotVariable(inc,indX, rangeA, num_concepts_mip, [.3 .3 .3] , condition, [0,4])
        function2use(i,indX, rangeA, num_concepts_mip, gcolor, condition, [0,4])
        ylabel('#concepts');
    subplot(nrows,ncols,6)
        plotVariable(inc, indX, rangeB, big_phi_mip, [.3 .3 .3] , condition, [0,0.6])
        function2use(i, indX, rangeB, big_phi_mip, gcolor, condition, [0,0.6])
        ylabel('$\Phi^{\rm Max}$','Interpreter','latex','FontSize',12);
        xlabel('#Generations');
end


%% I_SMMI and I_Pred (TODO)
if I_sen_mot == 1
    for i = inc
        figure(i+100)
        load(strcat(condition,'_ISenMot'));
        I_SMMI = results.ISenMot;
        HMot = results.HMot;
        HSen = results.HSen;

        load(strcat(condition,'_IPred'));
        I_Pred = results.IPred;
        H_States = results.HStates;

        figure(102)
        subplot(6,4,1+condNum)
            function2use(i, indX, rangeA, Fitness_level./(MaxFitness/100), gcolor, condition, [40,100])
            ylabel('Fitness (%)');
        subplot(6,4,5+condNum)
            function2use(i, indX, rangeA, I_SMMI, gcolor, condition, [0,1])
            ylabel('I_{SMMI}')
        subplot(6,4,9+condNum)
            function2use(i, indX, rangeA, HSen, gcolor, condition, [0,2])
            ylabel('H_{Sen}')   
        subplot(6,4,13+condNum)
            function2use(i, indX, rangeA, HMot, gcolor, condition, [0,2])
            ylabel('H_{Mot}')
        subplot(6,4,17+condNum)
            function2use(i, indX, rangeA, I_Pred, gcolor, condition, [0,3])
            ylabel('I_{Pred}')
        subplot(6,4,21+condNum)
            function2use(i, indX, rangeB, H_States, gcolor, condition, [0,4])
            ylabel('H_{State}')
            xlabel('#Generations');
    end
end

for i = inc
    if plot_network == 1
        g = rangeB(end);
        path = strcat(DPath, condition, '/trial');
        numNodes = 8;
        numMot = 2;
        numSen = 2;
        J_tempfile = strcat(path, int2str(evaluatedTrials(i)),'_', int2str(g), '_EdgeList.txt')
        J_temp = load(J_tempfile);

        if ~isempty(J_temp)
            J_temp = unique(J_temp, 'rows')+1;
            % MOTORS ARE SET TO 0 -> they don't actually have recurrent
            % connections
            J_temp = J_temp(J_temp(:,1) <= numNodes-numMot,:); 
            %Sensors might have incoming connection, but because they actually don't do anything they shouldn't be taken into account
            J_temp = J_temp(J_temp(:,2) > numSen,:);
            J_temp(J_temp(:,1) == J_temp(:,2),:)
            J_sparse = sparse(J_temp(:,1), J_temp(:,2),1,numNodes,numNodes);
            view(biograph(J_sparse))
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