function [Tx_Pos_total,Rx_Pos_total,Tx_Pos_total_free,Rx_Pos_total_free,temp_Tx_Pos_total,temp_Rx_Pos_total] = myTotalPosition(Tx_Pos,Rx_Pos,Tx_Pos_free,Rx_Pos_free)
% Preprocess to aggregate sensor positions according to the access order of Rx00-09 from Tx00-09 % 

for i = 1 : 24    
    temp_Tx = Tx_Pos(1 + (i-1)*10: i*10,:);
    temp_Tx_free = Tx_Pos_free(1 + (i-1)*10: i*10,:);
    temp_Rx = Rx_Pos(1 + (i-1)*10: i*10,:);
    temp_Rx_free = Rx_Pos_free(1 + (i-1)*10: i*10,:);
    for j = 1 : 10     
        temp_Tx_Pos_total(1 + (i-1)*100 + (j-1)*10:(i-1)*100 + j*10,:) = temp_Tx;
        temp_Rx_Pos_total(1 + (i-1)*100 + (j-1)*10:(i-1)*100 + j*10,:) = ones(10,1) * temp_Rx(j,:);
        temp_Tx_Pos_total_free(1 + (i-1)*100 + (j-1)*10:(i-1)*100 + j*10,:) = temp_Tx_free;
        temp_Rx_Pos_total_free(1 + (i-1)*100 + (j-1)*10:(i-1)*100 + j*10,:) = ones(10,1) * temp_Rx_free(j,:);
    end
end

% Remove the tx/rx positions from 580-600 (unobserved due to error)
Tx_Pos_total = [temp_Tx_Pos_total(1:580,:);temp_Tx_Pos_total(601:end,:)];
Rx_Pos_total = [temp_Rx_Pos_total(1:580,:);temp_Rx_Pos_total(601:end,:)];
Tx_Pos_total_free = [temp_Tx_Pos_total_free(1:580,:);temp_Tx_Pos_total_free(601:end,:)];
Rx_Pos_total_free = [temp_Rx_Pos_total_free(1:580,:);temp_Rx_Pos_total_free(601:end,:)];

end