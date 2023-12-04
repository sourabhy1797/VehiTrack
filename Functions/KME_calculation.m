function result = KME_calculation(x,y)
sigma = 1;
% disp(size(x));
% disp(size(y));
result = exp(-((x - y).^2) / (2 * sigma^2));
end
