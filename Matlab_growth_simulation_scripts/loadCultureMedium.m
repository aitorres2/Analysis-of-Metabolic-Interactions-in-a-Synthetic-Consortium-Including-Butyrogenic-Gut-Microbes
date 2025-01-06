% Created by Tomás Pröschle, 2022. Department of Chemical and Bioprocess Engineering, School of Engineering, Pontificia Universidad Católica de Chile, Santiago, Chile.
function [modelWithMedium, exchangesAdded] = ...
                    loadCultureMedium(model, cultureMedium, carbonSource)

% Load medium data
cultureMedium = readtable(cultureMedium);

% Set exchange reactions
[exchangeRxns, ~] = findExcRxns(model);
exchangeRxnsIDs = model.rxns(exchangeRxns);
model = changeRxnBounds(model, exchangeRxnsIDs, 0, 'l');
model = changeRxnBounds(model, cultureMedium.rxnID, ...
                              -cultureMedium.rates, 'l');

% Set carbon source
model = changeRxnBounds(model, carbonSource, -1, 'l');

% Check if there is growth
sol = optimizeCbModel(model, 'max');

thereIsGrowth = sol.f > 0;

if thereIsGrowth
    % Nothing more to add to the model
    exchangesAdded = {};
    modelWithMedium = model;

elseif ~thereIsGrowth
    % Look for other exchange reactions to add
    exchangeRxnsIDsNot = exchangeRxnsIDs;
    exchangeRxnsIDsNot(ismember(exchangeRxnsIDsNot, ...
                                cultureMedium.rxnID)) = [];
    exchangeCandidates = exchangeRxnsIDsNot;
    excSize = size(exchangeRxnsIDsNot, 1);
    div = excSize;
    
    while div > 0
        % Search by separating the set in halves
        div = fix(div/2);
        exchangeCandidatesInf = exchangeCandidates(1:div);
        exchangeCandidatesSup = exchangeCandidates(div + 1:end);

        % Inferior halve
        modelTest = changeRxnBounds(model, exchangeCandidatesInf, -1, 'l');
        solTestInf = optimizeCbModel(modelTest, 'max');

        % Superior halve
        modelTest = changeRxnBounds(model, exchangeCandidatesSup, -1, 'l');
        solTestSup = optimizeCbModel(modelTest, 'max');

        % Decide which halve to use
        if solTestInf.f > solTestSup.f
            exchangeCandidates = exchangeCandidatesInf;

        elseif solTestSup.f > solTestInf.f
            exchangeCandidates = exchangeCandidatesSup;

        elseif (solTestInf.f == 0) && (solTestSup.f == 0)
            break
        end
    end
    
    % With the new candidates, look for the minimum set that generates
    % growth    
    modelTestSearch = changeRxnBounds(model, exchangeCandidates, -1, 'l');
    solTestSearch = optimizeCbModel(modelTestSearch, 'max');

    if solTestSearch.f > 0
        % Create variables needed for the search
        searchForMinimumSet = true;
        excSize = size(exchangeCandidates, 1);
        excToTest = exchangeCandidates;
        excToKeep = {};
        excToRemove = {};
        keepCount = 1;
        removeCount = 1;

    else
        % If there isn't growth, then look for missing mets
        searchForMinimumSet = false;
        disp("Adding sink reactions to model")
        missingMets = biomassPrecursorCheck(modelTestSearch);

        % Initialize the model with the modelTestSearch
        modelWithMedium = modelTestSearch;

        % Iterate over the missingMets cell array
        for i = 1:length(missingMets)
            % Add a sink reaction for the current missing metabolite
            modelWithMedium = addSinkReactions(modelWithMedium, missingMets(i), -0.005, 1000); % 0.005 original - MODIFICAR AQUI PARA BOUND DE RXNS SINK
        end

        % Extract the newly added exchange reactions
        exchangesAdded = modelWithMedium.rxns(size(modelTestSearch.rxns, 1) + 1:end);
    end


    while searchForMinimumSet
        solList = zeros(excSize, 1);

        for i = 1:length(excToTest)
            % Remove the exchange reaction and test if it affects growth
            modelTest = changeRxnBounds(modelTestSearch, ...
                                                    excToTest{i}, 0, 'l');
            solTest = optimizeCbModel(modelTest);
            solList(i, 1) = solTest.f;

            if solTest.f > 0
                % If there is growth, remove
                excToRemove{removeCount, 1} = excToTest{i};
                removeCount = removeCount + 1;

            else
                % If there isn't growth, keep
                excToKeep{keepCount, 1} = excToTest{i};
                keepCount = keepCount + 1;
            end
        end
        
        % Pick new candidates
        exchangeCandidates = excToKeep;
        
        % Exit condition, all candidates are needed
        if all(solList == 0)
            % Candidates are selected
            exchangesAdded = exchangeCandidates;
            modelWithMedium = model;
            
            % Model is modified
            for i = 1:length(exchangesAdded)
                if startsWith(exchangesAdded{i}, 'sink_')
                    modelWithMedium = changeRxnBounds(modelWithMedium, ...
                                              exchangesAdded, -0.005, 'l'); % 0.005 original - MODIFICAR AQUI PARA BOUND DE RXNS SINK
                else
                    modelWithMedium = changeRxnBounds(modelWithMedium, ...
                                                  exchangesAdded, -1, 'l');
                end
            end
            break

        else
            % Keep iterating
            excSize = size(exchangeCandidates, 1);
            excToTest = exchangeCandidates;
            keepCount = 1;
        end
    end
end

end
