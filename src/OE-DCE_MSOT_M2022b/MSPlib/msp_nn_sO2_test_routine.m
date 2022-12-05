% load simulated data
load 'NNlib\sim_data_nn_001_AU151120'

% choose WLs
WL=[700:10:890];

% reshape to image stack
Recon = reshape(X0,1000,1000,20);
sO2_sim = reshape(y,1000,1000);

% unmix
sO2 = msp_nn_sO2(Recon,WL);

% compare
err= (sO2_sim-sO2);

% plot
figure
hist(err(:),100)



