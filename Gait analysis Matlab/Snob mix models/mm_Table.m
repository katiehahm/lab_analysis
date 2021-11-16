function mm_Table(mm)

%% How many classes and models per class
nClasses = mm.nClasses;
nModels  = mm.nModelTypes;

%% Print table header
fprintf('%10s', 'CLASS');
for k = 1:nClasses
    str = sprintf('%d (P=%.1f%%)', k, mm.a(k)*100);
    fprintf('%20s', str);
end
fprintf('%20s%20s', 'Model', 'Parameters');
fprintf('\n');

%% Print each model line
for i = 1:nModels    
        for k = 1:nClasses            
            model = mm.class{k}.model{i};           % model                
            theta = mm.class{k}.model{i}.theta;
            Ivar  = mm.class{k}.model{i}.Ivar;
            
            if(k == 1)
                str = mm.opts.VarNames{Ivar(1)};
                if(length(Ivar) > 1)
                    str = [str, '+'];
                end
                fprintf('%10s', str);                         
            end
            
            % Model types
            switch model.type                         
            case 'beta'
                str = sprintf('[%.2f %.2f]', theta(1),theta(2));                     
                mdlType = 'Beta';
            case 'Gaussian'    
                str = sprintf('[%.2f %.2f]', theta(1),sqrt(theta(2)));                
                mdlType = 'Gaussian';
            case 'multi'                  
                str = sprintf('[');
                for j = 1:length(theta)
                    str = [str, sprintf('%4.2f', theta(j))];
                    if(j ~= length(theta))
                        str = [str,' '];
                    end
                end
                str = [str,']'];             
                mdlType = 'Multinomial';                
            case 'weibull'
                str = sprintf('[%.2f %.2f]', theta(1),theta(2));                
                mdlType = 'Weibull';                
            case 'cfixweibull'
                str = sprintf('[%.2f %.2f]', theta(1),theta(2));                
                mdlType = 'Weibull Type I';
            case 'crndexp'
                str = sprintf('[%.2f %.2f]', theta(1),theta(2));                
                mdlType = 'Exp Rnd';
            case 'cfixexp'
                str = sprintf('%.2f', theta(1));                
                mdlType = 'Exp Type I';
            case 'exp'
                str = sprintf('%.2f', theta(1));
                mdlType = 'Exp';       
            case 'Laplace'
                str = sprintf('[%.2f %.2f]', theta(1),theta(2));                
                mdlType = 'Laplace';                
            case 'gamma'
                str = sprintf('[%.2f %.2f]', theta(1),theta(2));                        
                mdlType = 'Gamma';
            case 'Poisson'
                str = sprintf('%.2f', theta(1));                
                mdlType = 'Poisson';                
            case 'geometric'
                str = sprintf('%.2f', theta(1));                
                mdlType = 'Geometric';
            case 'invGaussian'
                str = sprintf('[%.2f %.2f]', theta(1),theta(2));                
                mdlType = 'Inv-Gaussian';  
            case 'negb'
                mu = theta(1); phi = theta(2);
                theta = [phi, 1-mu/(mu+phi)];
                str = sprintf('[%.2f %.2f]', theta(1),theta(2));                        
            	mdlType = 'Neg-binomial';                                            
            case 'vmf'
                str = '--'; 
                mdlType = 'von Mises-Fisher';        
            case 'linreg'
                str = '--'; 
                mdlType = 'Linreg';                     
            case 'logreg'
                str = '--'; 
                mdlType = 'Logreg';                                     
            case 'mvg'
                str = '--'; 
                mdlType = 'MVGaussian';                    
            case 'sfa'
                str = '--'; 
                mdlType = 'SingleFA';                                    
            otherwise
                   error('Unknown model type');
            end 
            fprintf('%20s', str);                      
        end
  
        switch model.type
            case 'beta'
                str = '[a, b]';
            case 'multi'
                str = '[proportions]';
            case 'Gaussian'
                str = '[mean sd]';
            case 'weibull'
                str = '[lambda k]';
            case 'cfixweibull'     
                str = '[lambda k]';
            case 'crndexp'                
                str = '[a,b]';
            case 'cfixexp'
                str = 'mean';
            case 'exp'
                str = 'mean';
            case 'Laplace'
                str = '[mu b]';
            case 'gamma'
                str = '[mu phi]';
            case 'Poisson'
                str = 'mean';
            case 'geometric'
                str = 'p';
            case 'invGaussian'                
                str = '[mu lambda]';
            case 'negb'
                str = '[r p]';                                
            otherwise
                str = '--';
        end
        fprintf('%20s%20s\n', mdlType, str);
    

end

end