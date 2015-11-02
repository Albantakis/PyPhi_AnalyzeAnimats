function Animat_SMMI
tic
clear all

Elem = 0:7;
numSen = 2;
Mot = [7 8];
%plotflag = 3;
trialnum = [0:399];

numtrials = numel(trialnum);
totsteps = 60000-1;
step = 512;
range = [0:step:totsteps];

cond = 'c3a2_change_c36a45';
%cond = 'c1a3_12sen_36';
%DPath = '~/Documents/Arend_XCodeAnimat2/temporalSpatialIntegrationLite/work_';
%DPath = '/Volumes/Macintosh HD 2/Simulations/Arend_XCodeAnimat/temporalSpatialIntegrationLite/work_';
DPath = '~/dev/animats/results/work_';
path = strcat(DPath, cond, '/trial');
% path3 = strcat(DPath2, cond,'_2', '/trial');
% path2 = strcat(DPath, cond, '_2', '/trial');
% cond2 = strcat(cond, '_2');
% cond3 = cond2;
path2 = path;
cond2 = cond;
%path2 = strcat(DPath2, cond, '_100','/trial');

count = 0;
evaluatedTrials = [];
for t = 1:numtrials
    Foldername = strcat(cond, '_trial', int2str(trialnum(t)));
    if exist(Foldername,'dir') == 7 
        count = count+1;
        evaluatedTrials = [evaluatedTrials trialnum(t)];
        if trialnum(t) > 99
            APath = strcat(path2, int2str(trialnum(t)), '_');
             if trialnum(t) > 599
                APath = strcat(path3, int2str(trialnum(t)), '_');
                cond = cond3;
             end
        else 
            APath = strcat(path, int2str(trialnum(t)), '_');
        end  
        for i = 1:length(range)
            %------------- get Fitness from Animat files---------------------------
             docname2 = strcat(APath, int2str(range(i)), '_KOdata.txt');
             Fitness = load(docname2);
             Fitness_level(count,i) = Fitness(1);  
            %Calculate probabilities for Sensory Motor predictive Information
            [p_Sen_t0, p_Mot_t1, p_Joint_t0t1, ~, ~] = Animat_LifeTimeDist(range(i), APath, numSen, Mot);
            p_Prod_t0t1 = p_Sen_t0*p_Mot_t1';   %This works if the ps are column vectors and states are ordered 000 100 010 110 etc
            p_Prod_t0t1 = p_Prod_t0t1(:);
            ISenMot(count,i) = KLD(p_Joint_t0t1, p_Prod_t0t1); %Mutual information of Sensor and Motors
            p_Mot_t1(p_Mot_t1 == 0) = 1;
            HMot(count,i) = DistEntropy(p_Mot_t1);
            p_Sen_t0(p_Sen_t0 == 0) = 1;
            HSen(count,i) = DistEntropy(p_Sen_t0);
            Hunpred(count,i) = HMot(count,i)-ISenMot(count,i);
        end
    end
end

figure
subplot(2,1,1)
 hold on
    plot(mean(Fitness_level), mean(ISenMot), '+r')
    xlim([1, 128])
subplot(2,1,2)
    hold on
    plot(range, mean(ISenMot), '-r')
    plot(range, mean(HMot), '-k')
    plot(range, mean(Hunpred), '-b')
    xlim([1, max(range)])

%------------------------ Save results ------------------------------------ 
results.Fitness = Fitness_level;
results.range = range;
results.ISenMot = ISenMot;
results.HMot = HMot;
results.HSen = HSen;
results.Hunpred = Hunpred;

save(strcat(cond, '_ISenMot'), 'results');
   
toc
end

function [p_InStates, p_OutStates, p_JointStates, GenpON, SizeLife] = Animat_LifeTimeDist(gen, APath, numSen, Mot)  
    fullstatefile = strcat(APath, int2str(gen), '_', 'Lifetime' ,'LogicTable.txt');
    Full_tpm = importdata(fullstatefile,',', 1);
    % Input Probability
    %--------------------------------------------
    Input_vec = Full_tpm.data(:,1:numSen);
    GenpON = sum(Input_vec)./size(Input_vec,1);
    [InStates, ~, inMap] = unique(Input_vec, 'rows', 'First');
    for i = 1:size(InStates, 1)
        %Order input States according to place in tpm
        InInd(i) = state2index(InStates(i,:), 2.*ones(size(InStates(1,:))));
    end
    inMap = InInd(inMap);
    InStates_distr = hist(inMap, 1:2^numSen)';
    p_InStates = InStates_distr./size(Input_vec,1);
    
    % Motor Probability
    %--------------------------------------------
    MotData = Full_tpm.data(:,Mot+9);
    SizeLife = size(MotData, 1);
    [OutStates, ~, outMap] = unique(MotData, 'rows', 'First');
    for i = 1:size(OutStates, 1)
        %Order output States according to place in tpm
        OutInd(i) = state2index(OutStates(i,:), 2.*ones(size(OutStates(1,:))));
    end
    outMap = OutInd(outMap);
      
    %Translate number of states into stateindex and make histogram
    OutStates_distr = hist(outMap, 1:2^size(MotData,2))';
    p_OutStates = OutStates_distr./sum(OutStates_distr,1);
    
    %Joint Probability
    %--------------------------------------------
    Joint_tpm = [Input_vec MotData]; %1ts
    %Joint_tpm = [Input_vec(1:end-1,:) MotData(2:end,:)];
    [JointStates, ~, JointMap] = unique(Joint_tpm, 'rows', 'First');
    for i = 1:size(JointStates, 1)
        %Order output States according to place in tpm
        JointInd(i) = state2index(JointStates(i,:), 2.*ones(size(JointStates(1,:))));
    end
    JointMap = JointInd(JointMap);
      
    %Translate number of states into stateindex and make histogram
    JointStates_distr = hist(JointMap, 1:2^size(Joint_tpm,2))';
    p_JointStates = JointStates_distr./sum(JointStates_distr,1);
end    

function entropy = DistEntropy(Dist)
    entropy = -sum(Dist.*log2(Dist));
end

function index = state2index(state_vec, num_states_vec)
index = state_vec(1) + 1;
for i = 2:length(num_states_vec)
    index = index + state_vec(i)*prod(num_states_vec(1:i-1));  
end
end