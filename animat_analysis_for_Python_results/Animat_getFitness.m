clear all
Elem = 0:7;
plotflag = 1;
%trialnum = [3 4 15 24 27 35 37 45 49 50 51 53 56 59 67 72 76 81 82 83 84 89 90 95 100 103 104 118 130 131 132 136 145 147 157 163 170 175 182 186 187 191 193 195 197];
%trialnum = [0 3 4 5 6 10 18 27 28 29 31 32 35 38 47 48 49 52 57 58 61 62 64 65 66 68 70 72 74 76 79 83 85 91 92 96 100 110 111 118 125 129 132 133 134 135 141 143 150 152 153 157 158 162 164 167 171 174 180 185 186 192 199];
%trialnum = [0 2 4 7 8 10 12 19 21 23 33 40 43 44 45 48 49 50 52 54 61 70 83 86 87 88 93 109 111 112 115 117 123 126 132 138 147 148 152 159 162 167 170 172 182 191 193 197];
%trialnum = [24 36 47 78 82 92 102 124 137 141 156 163 168 206 215 232 260 264 292 331 346 349 363 390]; %128 and 126
%trialnum = [24 36 47 92 102 137 141 156 163 215 232 264 292 331 346 349 363 390];

trialnum = 0:199;
numtrials = numel(trialnum);
totsteps = 30000-1;
step = 512;
range = 0:step:totsteps;
cond = 'c2a4_36';
%cond = 'c1a3_n001_rep20_36';
%cond = 'c1a3_n001_rep20_36'
%cond = 'c1a3_12Sen_36';
%cond = 'c1a3_n001_90000_36';
%cond = 'c1b11a2b5_36';
%cond = 'c23a45_36';
%cond = 'c35a271_36';
%cond = 'c1a3_change_c23a14';
%cond = 'c1a3_12sen_36';

%DPath = '~/Documents/Arend_XCodeAnimat2/temporalSpatialIntegrationLite/work_';
%DPath = '/Volumes/Macintosh HD 2/Simulations/Arend_XCodeAnimat/temporalSpatialIntegrationLite/work_';
DPath = '~/dev/animats/results/work_';
%DPath = '/Volumes/Macintosh HD 2/Simulations/AnimatsPaper2_data/animat_raw_files/work_';
path = strcat(DPath, cond, '/trial');
% path3 = strcat(DPath2, cond,'_2', '/trial');
%path2 = strcat(DPath2, cond, '/trial');
% cond2 = strcat(cond, '_2');
% cond3 = cond2;
path2 = path;
cond2 = cond;
%path2 = strcat(DPath2, cond, '_100','/trial');

MaxFitness = 128;

ElemUsed_temp = zeros(1, numel(Elem));
FitPhiCorr = [];

count = 0;
evaluatedTrials = [];
for t = 1:numtrials
   count = count+1;
   for i = 1:length(range)
        %------------- get Fitness from Animat files---------------------------
        if trialnum(t) > 99
            path = path2;
        end
        %docname = strcat(path, int2str(trialnum(t)), '_', int2str(range(i)), '_EdgeList.txt');
        docname2 = strcat(path, int2str(trialnum(t)), '_', int2str(range(i)), '_KOdata.txt');

        Fitness = load(docname2);
        Fitness_level(count,i) = Fitness(1);  
    end
end    
%cd ..
%%
File = strcat(cond, '_fitness');
save(File, 'Fitness_level', 'cond', 'range');

% if plotflag == 1
%     figure
%     subplot(3,1,1)
%     hold on
%     plot(range, MeanFitness)
%     xlim([1, totsteps])
% end