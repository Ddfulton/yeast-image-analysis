function [tau,MSD] = MSD_bead(dataname)

%% This code calculates the MSD versus lagtime (tau)  
 % input variable is dataname, for example 'dataout_2_1.mat'.
 % The quotation marks are important.
 
%% Example
 %  [tau,MSD] = MSD_bead( 'dataout_2_1.mat');
 %  to plot the MSD of all beads in chain 1:
 %  msd_1(:,:)=MSD(:,1,:);
 %  loglog(tau,msd_1)
 %  Note that the fixed beads are going to shift the axis a lot.  If you
 %  don't want to plot the MSd for bead 1, then
 %  clear msd_1
 %  msd_1(:,:)=MSD(2:52,1,:);
 %  loglog(tau,msd_1)
 %
 %  same thing if you just want to plot 1 bead (for example bead 26, chain 1)
 %  clear msd_1
 %  msd_1(:,:)=MSD(26,1,:);
 %  loglog(tau,msd_1)
 %
 %  All beads 26 from the 4 chains
 %  clear msd_1
 %  msd_1(:,:)=MSD(26,:,:);
 %  loglog(tau,msd_1)
 %  legend('Ch1', 'Ch2', 'Ch3', 'Ch4')
 
%% load data 
   load(dataname)
   
%% Subtract spb coordinates from each coordinate.
   xbc  = boundaryData.xbc;
   ybc  = boundaryData.ybc;
   Xspb = xbc(126,:);         Yspb = ybc(126,:);
   xsp  =  repmat(Xspb,52,1); ysp  = repmat(Yspb,52,1);
   xnew = 0*xstore;           ynew = xnew;
   for k=1:4
       tempx(:,:)  = xstore(:,k,:);  tempy(:,:)  = ystore(:,k,:);
       xnew(:,k,:) = tempx - xsp;    ynew(:,k,:) = tempy - ysp;
   end
   
%% Construct lagtime vector
   dt = t(2)-t(1);
   Tt = log(length(t))/log(10);
   temp = logspace(0, Tt, 100); % you can change the 100 if you want more
                                % points
   tau = unique(round(temp));   % tau is always going to have less entries 
                                % than temp      

%% Calculate MSD
   for k=1:4
       x(:,:) = xnew(:,k,:);   y(:,:) = ynew(:,k,:);
       for j=1:length(tau)
           dx  = x(:,1+tau(j):end) - x(:,1:end-tau(j));
           dy  = y(:,1+tau(j):end) - y(:,1:end-tau(j));
           dr2 = dx.^2 + dy.^2;
           MSD(:,k,j) = mean(dr2,2);     
       end
   end
   
   tau = tau*dt;