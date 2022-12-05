function r2 = rsquared(data,model)

m = mean(data);
SStot = sum((data-m).^2);
SSres = sum((data-model).^2);

r2 = 1 - SSres/SStot;