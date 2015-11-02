function Animat_HistogramsAtChange
% plot Fitness against each other 
changeGenerationNum = 59;
numBins = 14;
ISMMI_IPred = 0;
Hamming = 1;
TakeAverage = 1;
scatterhist = 0;

Task1 = {'c1a2', 'c1a3', 'c3a1'};
Task2 = {'c14a23'; 'c36a45'; 'c23a14'};

% Task1 = {'c3a1'};
% Task2 = {'c23a14'; 'c36a45'; 'c14a23'};

allCondCorrMatrix = [];
for j = 1:length(Task1)
    clear condition
    clear results
    if j == 3
        Task2 = {'c23a14'; 'c36a45'; 'c14a23'};
    end
    for i = 1:length(Task2)
        condition{i} = char(strcat(Task1{j}, '_change_', Task2{i}));
        results{i} = load(char(strcat(condition{i}, '_results')));
        
        if Hamming == 1
            cd ..
            cd(strcat('Hamming/', Task1{j}))

            HCM_prunedTrue{i} = load(strcat('HammingDistCM_', condition{i}, '_118_prunedTrue'));
            HTPM_prunedTrue{i} = load(strcat('HammingDistTPM_', condition{i}, '_118_prunedTrue'));

            cd ..
            cd ..
            cd('weighted')
        end
    end

    if ISMMI_IPred == 1
        %old!
        %Load only condA, cause all I am correlating with is the values at the
        %changeGenerationNum, which are the same for both conditions.
        load(strcat(condA, '_ISenMot'));
        ISMMIA = results;
        load(strcat(condA, '_IPred'));
        IPredA = results;
    end
    
    range1 = results{1}.range./10000;
    fitnessBeforeChange = max(results{1}.Fitness_level(:,55:59),[],2);
    indFitness = find(fitnessBeforeChange == 128 & results{1}.big_phi_mip(:,59) > 0);
    %indFitness = find(fitnessBeforeChange == 128);

    %% Correlations
    if Hamming == 1
        matrixForCorrelations = zeros(4+length(Task2)*4,length(indFitness));
    else
        matrixForCorrelations = zeros(4+length(Task2)*2,length(indFitness));
    end
%     P = mean(results{1}.big_phi_mip(indFitness,(changeGenerationNum-4):changeGenerationNum), 2);
%     C = mean(results{1}.num_concepts_mip(indFitness,(changeGenerationNum-4):changeGenerationNum), 2);
%     ZC = mean(results{1}.num_concepts_whole(indFitness,(changeGenerationNum-4):changeGenerationNum), 2);
    
    P = results{1}.big_phi_mip(indFitness,changeGenerationNum);
    C = results{1}.num_concepts_mip(indFitness,changeGenerationNum);
    ZC = results{1}.num_concepts_whole(indFitness,changeGenerationNum);
    NConn = double(results{1}.num_connections(indFitness,changeGenerationNum));

    matrixForCorrelations(1:4,:) = [P'; C'; ZC'; NConn'];

    for i = 1:length(Task2)   
       % rows 5-7 Fitness after Change
       matrixForCorrelations(4+(2*i-1),:) = (results{i}.Fitness_level(indFitness,changeGenerationNum)-results{i}.Fitness_level(indFitness,changeGenerationNum+1));
       % rows 8-11 Fitness at the end
       matrixForCorrelations(4+(2*i),:) = results{i}.Fitness_level(indFitness,end);
       % slope from dip to maximum
       if Hamming == 1
            matrixForCorrelations(10+(2*i-1),:) = HCM_prunedTrue{i}.HammingDist(indFitness, changeGenerationNum+1);
            matrixForCorrelations(10+(2*i),:) = HTPM_prunedTrue{i}.HammingDist(indFitness, changeGenerationNum+1);
       end
       
       if ISMMI_IPred == 1
            %old
            ISMMI = ISMMIA.ISenMot(indFitness, changeGenerationNum);
            IPred = IPredA.IPred(indFitness, changeGenerationNum);
       end
    end
    allCondCorrMatrix = [allCondCorrMatrix, matrixForCorrelations];
end
if TakeAverage == 1
    % very hacky
    allCondCorrMatrix_new = zeros(4+(size(allCondCorrMatrix, 1)-4)/length(Task2),size(allCondCorrMatrix,2));
    allCondCorrMatrix_new(1:4,:) = allCondCorrMatrix(1:4,:);
    for i = 5:size(allCondCorrMatrix_new, 1)
        allCondCorrMatrix_new(i,:) = mean(allCondCorrMatrix([i, i+2, i+4],:),1);
    end        
    allCondCorrMatrix = allCondCorrMatrix_new;
end
    
ind = 1:size(allCondCorrMatrix, 2);
%ind = find(allCondCorrMatrix()>0);
if ISMMI_IPred == 1
    [CorrMat, pCorr] = corr([FA60(ind), FAend(ind), FB60(ind), FBend(ind), P(ind), C(ind), ZC(ind) NConn(ind) ISMMI(ind) IPred(ind)],'type','Spearman');
else
    [CorrMat, pCorr] = corr(allCondCorrMatrix(:,ind)','type','Spearman');
end

figure; imagesc(CorrMat)
figure; imagesc(pCorr*100)
figure; imagesc(pCorr*100 <= 0.05)

CorrMat
pCorr

%% Scatter Plots and Histograms
if scatterhist == 1
for i = 1:length(Task2)
    for j = 1:4
        maxValue(j) = max(allCondCorrMatrix(j,:));
    end    
    
    %Phi
    xPhi = [0 0.0001:maxValue(1)/numBins:maxValue(1)+maxValue(1)/numBins]
    [nPhi, binPhi] = histc(allCondCorrMatrix(1,:), xPhi);  

    %Complex Concepts
    xConcepts = 0:maxValue(2)/numBins:maxValue(2)+maxValue(2)/numBins;
    [nConcepts, binConcepts] = histc(allCondCorrMatrix(2,:), xConcepts);

    %Zombie Concepts
    xZConcepts = 0:maxValue(3)/numBins:maxValue(3)+maxValue(3)/numBins;
    [nZConcepts, binZConcepts] = histc(allCondCorrMatrix(3,:), xZConcepts);  

    %Number of Connections Zombie
    xZConnections = 0:maxValue(4)/numBins:maxValue(4)+maxValue(4)/numBins;
    [nZConnections, binZConnections] = histc(allCondCorrMatrix(4,:), xZConnections);  

    if ISMMI_IPred == 1
        %I_SMMI
        maxISMMI = max(ISMMIA.ISenMot(indFitness, changeGenerationNum));
        xISMMI = 0:maxISMMI/numBins:maxISMMI+maxISMMI/numBins;
        [nISMMI, binISMMI] = histc(ISMMIA.ISenMot(indFitness, changeGenerationNum), xISMMI);

        %I_Pred
        maxIPred = max(IPredA.IPred(indFitness, changeGenerationNum));
        xIPred = 2:(maxIPred-2)/numBins:maxIPred+(maxIPred-2)/numBins;
        [nIPred, binIPred] = histc(IPredA.IPred(indFitness, changeGenerationNum), xIPred);
    end

figure
    Fitness60 = allCondCorrMatrix(4+(2*i-1),:);
    FitnessEnd = allCondCorrMatrix(4+(2*i),:);
    subplot(3,2,1)
    hold on
    plotFitness(xPhi, nPhi, binPhi, Fitness60, FitnessEnd, [1, 0, 0])
    xlim([-0.01 maxValue(1)+0.1])
    title('Phi')
    subplot(3,2,2)
    hold on
    plotFitness(xConcepts, nConcepts, binConcepts, Fitness60, FitnessEnd, [1, 0, 0])
    xlim([-0.5 maxValue(2)+0.5])
    title('#integrated concepts')
    subplot(3,2,3)
    hold on
    plotFitness(xZConcepts, nZConcepts, binZConcepts, Fitness60, FitnessEnd, [1, 0, 0])
    xlim([-0.5 maxValue(3)+0.5])
    title('#all concepts')
    subplot(3,2,4)
    hold on
    plotFitness(xZConnections, nZConnections, binZConnections, Fitness60, FitnessEnd, [1, 0, 0])
    xlim([-0.5 maxValue(4)+0.5])
    title('#all Connections')
    %Todo:
    if ISMMI_IPred == 1
        subplot(3,2,5)
        hold on
        plotFitness(xISMMI, nISMMI, binISMMI, ComplexA.Fitness_level,ComplexB.Fitness_level,indFitness)
        xlim([-0.1 maxISMMI+0.1])
        title('I_{SMMI}')
        subplot(3,2,6)
        hold on
        plotFitness(xIPred, nIPred, binIPred, ComplexA.Fitness_level,ComplexB.Fitness_level,indFitness)
        xlim([1.9 maxIPred+0.1])
        title('I_{Pred}')
    end
end
end
end

function plotFitness(x, n, bin, Fitness60, FitnessEnd, color)
    bar(x, n)
    scatter(x(bin), Fitness60, 'o', 'MarkerEdgecolor', [0.5 0.5 0.5])
    scatter(x(bin), FitnessEnd, '*','MarkerEdgecolor', [0.5 0.5 0.5])
    [yA, errorA, indA] = meanBinValues(length(n),bin, Fitness60);
    errorbar(x(indA), yA(indA), errorA(indA), '+','color', color)
    [yB, errorB, indB] = meanBinValues(length(n),bin, FitnessEnd);
    errorbar(x(indB), yB(indB), errorB(indB), '+','color', [1 1 0])
end

function [result, stdResult, indnan] = meanBinValues(maxbin, bins, data)
    result = nan(maxbin,1);
    stdResult = nan(maxbin,1);
    for i = 1:maxbin
        ind = bins == i;
        result(i) = nanmean(data(ind));
        stdResult(i) = nanstd(data(ind));
    end
    indnan = find(~isnan(result));
end