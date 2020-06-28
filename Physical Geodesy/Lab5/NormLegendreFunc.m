function [P] = NormLegendreFunc(l, t)

n = length(t);
P = zeros(l+1, l+1, n);
P(1,1,:) = ones(1, n);
P(2,1,:) = sqrt(3)*t;
P(2,2,:) = sqrt(3*(1-t.^2));

for i = 2:1:l
    for j = 0:1:l
        if i == j
            P(i+1,j+1,:) = sqrt((2*i+1)/(2*i))*(P(2,2,:)/sqrt(3)).*P(i,j,:);
        else
            P(i+1,j+1,:) = (sqrt(2*i-1)*(P(2,1,:)/sqrt(3)).*P(i,j+1,:)-sqrt((i-1+j)*(i-1-j)/(2*i-3))*P(i-1,j+1,:))*sqrt((2*i+1)/(i-j)/(i+j));
        end
    end
end

end

