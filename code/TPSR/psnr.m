function decibels = psnr(A,B)
maxval = max(max(max(A)), max(max(B)));
minval = min(min(min(A)), min(min(B))); 
A = (A-minval)/(maxval-minval)*8;
B = (B-minval)/(maxval-minval)*8; 
error_diff = A - B;
decibels = 20*log10(1/(sqrt(mean(mean(error_diff.^2))))); 