function [J] = replicate(I)
    [M, N] = size(I);

    J = zeros(M+2,N+2);

    J(1,2:end-1) = I(1,:);
    J(end,2:end-1) = I(end,:);
    J(2:end-1,1) = I(:,1);
    J(2:end-1,end) = I(:,end);

    J(2:end-1, 2:end-1) = I;
end
