% plot Fitness against each other 
plotflag =3;
DPath = '/Users/larissa_alb/dev/animats/results/work_';
numNodes = 8;
numMot = 2; 
numSen = 2;

condi = 'c1a3'

condA = strcat(condi, '_change_c14a23');
condB = strcat(condi, '_change_c36a45');
pathA = strcat(DPath, condA, '/trial');
pathB = strcat(DPath, condB, '/trial');

load(strcat(condA, '_dataCB'));

index = evaluatedTrials;

range1 = range./10000;
range1S = 1:length(range);
%Zombie = load(strcat(condi, '_36_ZombiedataAllC'));
PhiMip1 = BigPhiMip;
NConceps = MeanNumConcepts;
Fitness1 = Fitness_level;
clear Fitness_level

load(strcat(condB, '_dataCB'));
PhiMip2 = BigPhiMip;
Fitness2 = Fitness_level;

if plotflag == 1
    % plot Fitness against each other
    figure
    hold on
    indP = [];
    indP0 = [];
    for i = 1:size(Fitness1,1)
        if max(BigPhiMip(i,55:59)) > 0
            indP = [indP i];
            plot(Fitness1(i,60:end)', Fitness2(i,60:end)', 'k')
            plot(Fitness1(i,60)', Fitness2(i,60)', 'ok')
            plot(Fitness1(i,end)', Fitness2(i,end)', '*k')
        else
            indP0 = [indP0 i];
            plot(Fitness1(i,60:end)', Fitness2(i,60:end)', 'b')
            plot(Fitness1(i,60)', Fitness2(i,60)', 'ob')
            plot(Fitness1(i,end)', Fitness2(i,end)', '*b')
        end
    end
    plot(mean(Fitness1(indP,60:end)), mean(Fitness2(indP,60:end)),'r')
    plot(mean(Fitness1(indP0,60:end)), mean(Fitness2(indP0,60:end)),'g')

elseif plotflag == 2
%% plot PhiMip against each other
    figure
    hold on
    indP = [];
    indP0 = [];
    for i = 1:size(PhiMip1,1)
        if max(BigPhiMip(i,55:59)) > 0
            indP = [indP i];
            plot(PhiMip1(i,60:end)', PhiMip2(i,60:end)', 'k')
            plot(PhiMip1(i,60)', PhiMip2(i,60)', 'ok')
            plot(PhiMip1(i,end)', PhiMip2(i,end)', '*k')
        else
            indP0 = [indP0 i];
            plot(PhiMip1(i,60:end)', PhiMip2(i,60:end)', 'b')
            plot(PhiMip1(i,60)', PhiMip2(i,60)', 'ob')
            plot(PhiMip1(i,end)', PhiMip2(i,end)', '*b')
        end
    end
    plot(mean(PhiMip1(indP,60:end)), mean(PhiMip2(indP,60:end)),'r')
    plot(mean(PhiMip1(indP0,60:end)), mean(PhiMip2(indP0,60:end)),'g')

elseif plotflag == 3
    %color = [max(PhiMip1(:,55:59),[],2), 0.5.*ones(size(Fitness1, 1),1), max(max(PhiMip1(:,55:59),[],2)) - max(PhiMip1(:,55:59),[],2)]./max(max(PhiMip1(:,55:59),[],2));
    color = [max(NConceps(:,55:59),[],2), 0.5.*ones(size(Fitness1, 1),1), max(max(NConceps(:,55:59),[],2)) - max(NConceps(:,55:59),[],2)]./max(max(NConceps(:,55:59),[],2));

    figure
    subplot(1,2,1)
    hold on
    for i = 1:size(Fitness1,1)
        if max(PhiMip1(i,55:59)) > 0
            plot(range, Fitness1(i,:), 'color', [color(i,:)])
        else
            plot(range, Fitness1(i,:), '-k')
        end
    end
    
    subplot(1,2,2)
    hold on
    for i = 1:size(Fitness2,1)
        if max(PhiMip1(i,55:59)) > 0
            plot(range, Fitness2(i,:), 'color', [color(i,:)])
        else
            plot(range, Fitness2(i,:), '-k')
        end
    end
elseif plotflag == 4
    totsteps = 60000-1;
    range = (0:512:totsteps);
    gen = range([59:60 end]);
    color = [max(NConceps(:,55:59),[],2), 0.5.*ones(size(Fitness1, 1),1), max(max(NConceps(:,55:59),[],2)) - max(NConceps(:,55:59),[],2)]./max(max(NConceps(:,55:59),[],2));
    
    trialnum = 35; %1:size(Fitness1,1);
    
    for t = trialnum
        figure(100+t)
        hold on
        plot(range, Fitness1(t,:), 'color', [color(t,:)])
        plot(range, Fitness2(t,:), 'color', [color(t,:)])
    
    
        for g = gen;
        %------------- get Fitness from Animat files---------------------------
        J_tempfile = strcat(pathA, int2str(index(t)),'_', int2str(g), '_EdgeList.txt');
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
        
        for g = gen(2:end);
        %------------- get Fitness from Animat files---------------------------
        J_tempfile = strcat(pathB, int2str(index(t)),'_', int2str(g), '_EdgeList.txt');
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