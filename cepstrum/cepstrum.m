% Your main loop
clear; clc; close all
addpath("voice-icar-federico-ii-database-1.0.0\")
numFiles = 208; % Total number of files
classes = {}; % To store all unique classes
quefrenciesByClass = struct();
top_k_Quefrencies = 7; 

for k = 1:numFiles
    % Generate filename for each iteration
    infoFilename = sprintf('voice%03d-info.txt', k);
    
    % Extract diagnosis and gender from the file
    info = importdata(infoFilename);
    diagnosis = info{5};
    diagnosis = diagnosis(12:end); % extract only the diagnosis part
    
    % Renaming diagnosis based on the conditions
    if contains(diagnosis, 'hyperkinetic dysphonia')
        diagnosis = 'hyperkinetic dysphonia';
    elseif contains(diagnosis, 'healthy')
        diagnosis = 'healthy';
    elseif contains(diagnosis, 'hypokinetic dysphonia')
        diagnosis = 'hypokinetic dysphonia';
    elseif contains(diagnosis, 'reflux laryngitis')
        diagnosis = 'reflux laryngitis';
    end
    diagnosisGender = strcat(diagnosis, '_', gender);
    
    % Add the diagnosis to the classes array if it's not already there
    if ~ismember(diagnosis, classes)
        classes{end+1} = diagnosis;
        quefrenciesByClass.(matlab.lang.makeValidName(diagnosis)) = [];
        topQuefrenciesByClass.(matlab.lang.makeValidName(diagnosis)) = [];
    end
    
    % Compute quefrencies for the voice 
    [voiceQuefrencies, shifts] = computeQuefrenciesForVoice(k);
    % Append the quefrencies to the corresponding class
    quefrenciesByClass.(matlab.lang.makeValidName(diagnosis)) = [quefrenciesByClass.(matlab.lang.makeValidName(diagnosis)) voiceQuefrencies];

    % Extract the top 6 quefrencies
    [voiceQuefrenciesTop, ~] = computeQuefrenciesForVoice(k, top_k_Quefrencies);
    topQuefrenciesByClass.(matlab.lang.makeValidName(diagnosis)) = [topQuefrenciesByClass.(matlab.lang.makeValidName(diagnosis)); voiceQuefrenciesTop];
end


% Compute and plot the mean and variance for each class
means = zeros(length(classes), length(voiceQuefrencies));
variances = zeros(length(classes), length(voiceQuefrencies));

folders = {'mean', 'variance', 'plots', 'histograms'};
for idx = 1:length(folders)
    if ~exist(folders{idx}, 'dir')
        mkdir(folders{idx});
    end
end

for i = 1:length(classes)
    classData = quefrenciesByClass.(matlab.lang.makeValidName(classes{i}));
    means(i, :) = mean(classData, 2);
    variances(i, :) = var(classData');

    figure('Visible', 'off');
    plot(shifts, means(i, :));
    title(['Mean Quefrencies for ', classes{i}]);
    ylabel('Mean');
    xlabel('Quefrency');
    grid on
    axis tight;

    % Save the mean figure in the 'mean' folder
    saveas(gcf, fullfile( 'mean', [classes{i}, '_mean.png']));

    % Plotting the variance for the class
    figure('Visible', 'off');
    plot(shifts, variances(i, :));
    title(['Variance of Quefrencies for ', classes{i}]);
    ylabel('Variance');
    xlabel('Quefrency');
    grid on
    axis tight;

    % Save the variance figure in the 'variance' folder
    saveas(gcf, fullfile('variance', [classes{i}, '_variance.png']));

end

% Generate joint histograms for each class and quefrency
for i = 1:top_k_Quefrencies - 2
    q1_m = topQuefrenciesByClass.healthy(i, :);
    q2_m = topQuefrenciesByClass.hyperkineticDysphonia(i, :);
    q3_m = topQuefrenciesByClass.hypokineticDysphonia(i, :);
    q4_m = topQuefrenciesByClass.refluxLaryngitis(i, :);

    figure('Visible', 'off'); % Create a new figure for this iteration

    subplot(2,2,1);
    histogram(abs(q1_m), 12, 'Normalization', 'probability', 'FaceAlpha', 0.5, 'DisplayName', 'Healthy - Male');
    title('Healthy');

    subplot(2,2,2);
    histogram(abs(q2_m), 12, 'Normalization', 'probability', 'FaceAlpha', 0.5, 'DisplayName', 'Hyperkinetic Dysphonia - Male');
    title('Hyperkinetic Dysphonia');

    subplot(2,2,3);
    histogram(abs(q3_m), 12, 'Normalization', 'probability', 'FaceAlpha', 0.5, 'DisplayName', 'Hypokinetic Dysphonia - Male');
    title('Hypokinetic Dysphonia');

    subplot(2,2,4);
    histogram(abs(q4_m), 12, 'Normalization', 'probability', 'FaceAlpha', 0.5, 'DisplayName', 'Reflux Laryngitis - Male');
    title('Reflux Laryngitis');

    legend();

    sgtitle(['Quefrency ' num2str(i)]);
    saveas(gcf, fullfile('histograms', sprintf('Quefrency_m_%d.png', i)));

    
end
figure('Visible', 'off');
for i = 1:length(classes)
    subplot(4,1,i)
    classData = quefrenciesByClass.(matlab.lang.makeValidName(classes{i}));
    k = randi(size(classData, 2));
    plot(shifts, classData(:, k));
    title(['Cepstrum of ', classes{i}, ' Sample ', num2str(k) ]); 
    ylabel('C_{r}');
    xlabel('Quefrency');
    grid on
    axis tight;
    ylim([-1, 1])
end

saveas(gcf, fullfile(['plots', 'allclasses.png']));