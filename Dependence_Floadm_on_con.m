                    %% F_load_m  = function(con)
% 23.6.2025
% Ondočko Štefan - stefan.ondocko@tuke.sk

% SIMULINK MODEL Setting Examle
% 1) Open "Model_x_x"
% 2) Check the configuration of the Repeating Sequence Stair generator
%    -Vector of output values: [1e2 3e2 1e3 3e3 1e4 3e4 1e5 3e5 1e6 3e6 ].'
%    ..have 10 elements.
%    -Sample time: 10sec
% 4) Stop Time = "Number of elements in the output vector of the Repeating Sequence Stair block" x "Sample time"
%    Stop Time = 10 x 10sec = 100sec.
% 5) Run simulation of the "Model_x_x".
% 6) Run script "Dependence_Floadm_on_con.m"

% If the optimal constant "con" for the calculation is found (xx),
% the script "Dependence_Floadm_on_con.m" is no longer needed.
% The optimal calculation is then performed as follows:
% 7) Manual Switch to selected optimal Constant-"con" (xx)
% 8) Stop Time = 1 x 10sec = 10sec.

clc, close all
format long

% Check inputs
if exist('F_load_out','var') &&...
        exist('con_out','var') &&...
        exist('F_load_m_out','var') &&...
        exist('t_out','var') && sum(con_out)/size(con_out,1)~=con_out(1)

% >---------------------------- Inputs ------------------------------------
    %con_out "Decimation=1" in Simulink models!
CON=[con_out(1:size(con_out,1)-1)]';

    %Load force vector - "F_load"
F_load = F_load_out;

    %Control measurement of load force vector - "F_load_measured",all sequences
F_load_m = [F_load_m_out'; t_out'];

    %Lag, for the steady state
Lag=1.2;
    Lag_corr = Lag+(1e-12); %Not from the absolute zero

    %Lead, in the steady state
Lead=0.1;
% ------------------------------------------------------------------------<

SampleTime = round(max(t_out)/size(CON,2),0); %"100/10"
Sequence = size(CON,2); %"10s"
StopTime = SampleTime*Sequence; %"100s"

    %Sequences vector
SequenceVector = SampleTime*(0:1:StopTime/SampleTime);
    
    %Samples
plot(t_out,F_load_m_out,'x')
hold on 

for k=1:size(SequenceVector,2) %Cycle of the sequence
    if k < size(SequenceVector,2) %Without the last value of the vector "Time"
        for i=1:size(F_load_m,2) %Cycle of the selection
            if ((SequenceVector(k)+Lag_corr)<=F_load_m(2,i)) && (F_load_m(2,i)<=(SequenceVector(k+1)-Lead))
                F_load_m_separ(:,i) = F_load_m(:,i);
            end
        end %End of the selection cycle
        F_load_m_separ_k = F_load_m_separ;

        clear F_load_m_separ %Clear previous values
        F_load_m_separ_k; %Just for the check..

    F_load_m_separ_1 = F_load_m_separ_k(1,:); %Force for the sequence
    F_load_m_separ_2 = F_load_m_separ_k(2,:); %Time for the sequence

        %Combining data into Cells
        A{k}=F_load_m_separ_1(F_load_m_separ_1~=0); %Cell of k-th arrays without zeros
        n(k) = size(F_load_m_separ_1(F_load_m_separ_1~=0),2); %Number of the elements without zeros
        Max_F_load_m(k) = max(A{k}); %Maxims
            MAX_F_load_m{k} = Max_F_load_m(k).*ones(1,n(k)); %The array of the maxims for plotting
        Min_F_load_m(k) = min(A{k}); %Minims
            MIN_F_load_m{k} = Min_F_load_m(k).*ones(1,n(k)); %The array of the minims for plotting
    
        %Mean
    Mean_F_load_m(k) = sum(F_load_m_separ_1)./n(k); %Mean
        MEAN_F_load_m{k} = Mean_F_load_m(k).*ones(1,n(k)); %The array of the Mean

    % >------------------------- No necessery -----------------------------
        %Deviation from the arithmetic mean
    Delta_F_load_m = F_load_m_separ_1(F_load_m_separ_1~=0)-Mean_F_load_m(k);      
        
        if n(k)<=30
                %Standard error of the mean (n<=30)
            Sigma_med_F_load_m(:,k) = sqrt(sum(Delta_F_load_m.^2)/(n(k)*(n(k)-1)));

                %Standard deviation (n<=30)
            Sigma_F_load_m(:,k) = sqrt(sum(Delta_F_load_m.^2)/(n(k)-1));

            %SIGMA_F_load_m{k} = Sigma_F_load_m(k).*ones(1,n(k)); %The array of Sigmas for plotting
        else
                %Standard error of the mean
            Sigma_med_F_load_m(:,k) = sqrt(sum(Delta_F_load_m.^2)/n(k)^2);

                %Standard deviation
            Sigma_F_load_m(:,k) = sqrt(sum(Delta_F_load_m.^2)/n(k));

            %SIGMA_F_load_m{k} = Sigma_F_load_m(k).*ones(1,n(k)); %The array of Sigmas for plotting
        end
    % --------------------------------------------------------------------<

    time = F_load_m_separ_2(F_load_m_separ_2~=0); %Time without zeros
        TIME{k} = time; %The array of the times'
  
        R=plot(TIME{k}, MAX_F_load_m{k},"Color",'r','LineStyle','-','LineWidth',2);
        G=plot(TIME{k},MEAN_F_load_m{k},"Color",'g','LineStyle','-','LineWidth',3);
        BL=plot(TIME{k}, MIN_F_load_m{k},"Color",'bl','LineStyle','-','LineWidth',2);
        
        % >------------ No necessery. 3x Sigma for plotting ---------------
        % plot(TIME{k},MEAN_F_load_m{k}+3*SIGMA_F_load_m{k}, ...
        %      TIME{k},MEAN_F_load_m{k}-3*SIGMA_F_load_m{k},"Color",'g','LineStyle','--','LineWidth',2)
        % ----------------------------------------------------------------<
    else
    end

end %End of the sequences cycle

    %Data F_load
plot(t_out,F_load,"Color",'k','Marker','.','LineWidth',2)
hold off

    %Determining the limits of the graph
korekt=(min(Min_F_load_m)-max(Max_F_load_m))/10;
    if (min(Min_F_load_m) <= F_load) && (F_load <= max(Max_F_load_m))
ylim([min(Min_F_load_m) max(Max_F_load_m)-korekt]);
    elseif max(Max_F_load_m) < F_load
ylim([min(Min_F_load_m)+korekt F_load-korekt]);
korekt=10*korekt;
    else
ylim([F_load max(Max_F_load_m)]);
    end
xlim([0 StopTime]),grid
xticks(0:SampleTime:StopTime); %Custom-made grid according to x

        %Graph text of the value "con"
    Xtext=SequenceVector(1:(size(SequenceVector,2)-1))+1;
    Ytext=max(Max_F_load_m)-korekt/2;
    for h=1:size(Xtext,2)
        text(Xtext(h), Ytext, 0,{num2str(CON(h),'%.e')});
    end

    %Axis marking
title('Accuracy of "F_{load m}" calculation depending on "con"')
xlabel('t [s]')
ylabel('F_{load m} [N]')

    %Legend
dash=yline(0,'--',"Color",'k','LineStyle','-','LineWidth',2);
legend([dash R(1) G(1) BL(1)],...
    {'F_{load}','Max F_{load m}','Mean F_{load m}','Min F_{load m}'},...
    'Location','south');

n %List of the samples per sequnce
Mean_F_load_m %List of the Mean_F_load_m per sequnce
Sigma_F_load_m; %List of Sigma values if need it

elseif exist('con_out','var') &&...
        sum(con_out)/size(con_out,1)==con_out(1)
    disp('Probably "Manual Switch" in SIMULINK model is in wrong position')
else
    disp('Incomplete input DATA.. Run the Simulik Model')
end

% Ing. Ondočko Štefan, PhD.
% -------------------------------------------------------
% Katedra výrobnej techniky a robotiky
% Strojnícka fakulta
% Technická univerzita Košice
% Slovensko
% Telefón: 055/602 3238
% -------------------------------------------------------
% Department of Manufacturing Machinery and Robotics
% Faculty of Mechanical Engineering
% Technical University of Košice
% Slovakia
% Telephone: 055/602 3238
% -------------------------------------------------------
% www.sjf.tuke.sk
% https://orcid.org/0000-0002-2236-8792
% Scopus Author ID: 57212252729
