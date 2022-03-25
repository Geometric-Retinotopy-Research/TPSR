function R2 = calcR2(yhat, y)

SSE = sum((y-yhat).^2);

SST = sum((y - mean(y)).^2);

R2   = (1- SSE/SST)*100;


