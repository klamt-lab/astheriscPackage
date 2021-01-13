function [reportPath] = getPublicationRunReportPath(minimalYieldFactor, numMaxExchanges, modelPath)
    if contains(modelPath, "ecgsDouble")
        modelName = "ecgsDouble";
    elseif contains(modelPath, "iML1515double")
        modelName = "iML1515double";
    elseif contains(modelPath, "ecolicore2double")
        modelName = "ecolicore2double";
    elseif contains(modelPath, "ecolicore2triple")
        modelName = "ecolicore2triple";
    else
        error("No model name found :O")
    end

    reportPath = "./0Astherisc/publication_runs/ecoli_models/run_results/report__model_"+modelName+"__maxYield_"+minimalYieldFactor+"__numMaxExchanges_"+numMaxExchanges+".txt";
end
