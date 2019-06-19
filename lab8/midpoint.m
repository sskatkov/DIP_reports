function [image, snrs] = midpoint(I, J)
    [M, N] = size(I);
    image = zeros(M,N);
    for i = 2:M+1
        for j = 2:N+1
            block = J((i-1):(i+1), (j-1):(j+1));
            block = block(:);
            image(i-1,j-1) = (max(block) + min(block))/2;
        end
    end
    snrs = snr(double(I),double(image)-double(I));
end