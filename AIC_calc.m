function aic=AIC_calc(data)

    % recursive check for matrices
    if ~isvector(data)
        aic = zeros(size(data));
        for i = 1:size(data,2)
            aic(:,i) = AIC_calc(data(:,i));
        end
        return
    end
    %%
    N = length(data);
    aic = zeros(size(data));

%     for i = 1:length(data)
%        var1 = var(data(1:i)) ;
%        var2 = var(data(i+1:end));
% 
%        aic(i) = i*log10(var1) + (N-i-1)*log10(var2);
%        allVar2(i) = var2;
%        allVar1(i) = var1;
%        
% 
%     end
% 
    [varForward meanForward]=cumVar(data);
    varBack = cumVar(flipud(fliplr(data(2:end))));
    varBack = flipud(fliplr(varBack));
    varBack(end+1) = varBack(end);
   
    % this should be equivlent to aic but its not.... will speed up code.
    aic = log10(varForward).*[1:N]' + log10(varBack).*[N-2:-1:-1]';

    % remove the mean...
    aic  = aic- mean(aic(~isnan(aic)));
%     clf
%     subplot(3,1,1)
%     plot(allVar1)
%     hold on
%     plot(varForward)
%     
%     hold off
%     subplot(3,1,2)
%     plot(allVar2)
%     hold on
%     plot(varBack)
%     
%     subplot(3,1,3)
%     plot(aic)
%     hold on
%     plot(aic2)
%     hold off
    
end
%%
function outMean = cumMean(data)

Nm = [1:length(data)]';
outMean=cumsum(data)./Nm;

end
%%
function [outVar outMean] = cumVar(data)
    Nm = [1:length(data)]';
    outMean = cumMean(data);
    %outVar = abs(data - outMean).^2;
    %outVar = cumsum(outVar);
    %outVar  = outVar./(Nm-1);
    outVar  = cumsum(data.^2) - Nm.*(outMean).^2;
    outVar = outVar./(Nm-1);
end

%for nn=imin:(imax-1)
%         tic

%    AIC(nn)=nn*log(var(data(imin-500:nn)))+(N-nn-1)*log(var(data(nn+1:N)));

%end

% figure; plot(AIC)
%[val,I]=min(AIC);
