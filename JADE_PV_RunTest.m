clear all
clc

format long;
format compact;
warning('off');

runtime=50;

tic()
    indexModel=input('Input the type of PV model:');
    
    if indexModel==1
        maxNFEs=10000;
    elseif indexModel==2
        maxNFEs=20000;
    elseif indexModel==3
        maxNFEs=30000;
    elseif indexModel==4
        maxNFEs=30000;
    elseif indexModel==5
        maxNFEs=30000;
    end

    NP=50;
    for i=1:runtime
        [recordRMSE(i,:),bestP(i,:),bestFitness(i)]=JADE_PV_demo(indexModel,maxNFEs,NP);
    end
    
    ResultStatic=[std(bestFitness),max(bestFitness),min(bestFitness),mean(bestFitness)];
    [minRMSE,index]=min(bestFitness);
    ParameteResult=bestP(index,:);  
    plotRMSE=recordRMSE(index,:);
    
    disp(strcat('第',num2str(indexModel),'种模型结果为：'));
    
    if indexModel==1||indexModel==3||indexModel==4||indexModel==5
        disp(strcat('I_ph=',num2str(bestP(index,1)),',I_sd=',num2str(bestP(index,2)),',R_s=',num2str(bestP(index,3)),...
            ',R_sh=',num2str(bestP(index,4)),',a=',num2str(bestP(index,5)),',RMSE=',num2str(minRMSE)));
    elseif indexModel==2
        disp(strcat('I_ph=',num2str(bestP(index,1)),',I_sd1=',num2str(bestP(index,2)),',R_s=',num2str(bestP(index,3)),...
            ',R_sh=',num2str(bestP(index,4)),',a1=',num2str(bestP(index,5)),',I_sd2=',num2str(bestP(index,6)),...
            ',a2=',num2str(bestP(index,7)),',RMSE=',num2str(minRMSE)));
    end
    
    disp(strcat('Std=',num2str(ResultStatic(1)),',Max=',num2str(ResultStatic(2)),',Min=',num2str(ResultStatic(3)),...
        ',Mean=',num2str(ResultStatic(4))));
    
toc()

