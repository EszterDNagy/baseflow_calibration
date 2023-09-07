function [beta_0,Qb_0,Tc_0,Te_0,Tl_0,f] = beta_cal(t,Q,P,tpeak,A,ts,plot_mode)
%This function calibrates the parameter (beta) of the Lyne-Hollick (1979) recursive 
%digital filter for base flow separation.

%   INPUTS
%   t - datetime array of the discharge and precipitation time series
%   Q and P - Discharge and precipitation time series with uniform length and
%             timestep.
%   tpeak - Datetime array of runoff peaks. These identify the
%           rainfall-runoff events which will be used for calibration.
%   A - catchment area in km2
%   ts - timestep of the time series in seconds

%   OUTPUTS
%   beta_0 - Calibrated coefficient.
%   Qb_0 - Base flow time series calculated using beta_0.
%   Tc_0, Te_0, Tl_0 - time of concentration, time to equilibrium, and lag
%                      time using beta_0 for the calibration events.
%   f - Figure containing the plots of the main results.

% array of beta values to be tested
beta = [.925:.005:.99, .991:0.001:.999, .9991:.0001:.9999]';

% create empty arrays
t_qd_s = NaT(numel(beta),numel(tpeak));
t_qd_e = t_qd_s; t_p_e = t_qd_s; t_p_s = t_qd_s; t_P_com = t_qd_s; t_qd_com = t_qd_s;
Tc = zeros(numel(beta),numel(tpeak)); Te = Tc; Tl = Tc; r = Tc; FI = Tc;
P_com = Tc; qd_com = P_com;

% calculate time parameters (Tc, Te, Tl) for each beta and each event
for i = 1:numel(beta)
    % 1. base flow separation using the recursive filter
    Qb = qsep(Q,beta(i));
    Qd = Q - Qb;
    for j = 1:numel(tpeak)
        % 2. find start and end of runoff event using base flow
        % ID of peak discharge for given event
        id_peak = find(abs(hours(t-tpeak(j))) == min(abs(hours(t-tpeak(j)))),1,'first');%id_peak = find(t==tpeak(j));
        % find start and end point of runoff when base flow becomes larger
        % than 90% of the total flow befor and after the peak
        try
            t_qd_s(i,j) = t(find(Qb(1:id_peak) > .9.*Q(1:id_peak),1,'last')); 
        catch
            t_qd_s(i,j) = t(1);
        end
        try
            t_qd_e(i,j) = t(id_peak-1+find(Qb(id_peak:end) > .9.*Q(id_peak:end),1,'first'));
        catch
            t_qd_e(i,j) = t(end);
        end
        
        % 3. calculate total direct runoff using catchment area and timestep
        R = sum(Qd(find(t==t_qd_s(i,j)):find(t==t_qd_e(i,j))))/A*ts/1000;
        
        % 4. find start and end of precipitation
        % end of precipitation is 
        t_p_e(i,j) = t(find(P(1:id_peak)~=0,1,'last'));
        if t_p_e(i,j) == t(id_peak)
            t_p_e(i,j) = t(id_peak+find(P(id_peak:end)==0,1,'first')-2);
        end
        k = 0; P_event = 0; Pk = P(t==t_qd_s(i,j));
        while R > sum(P_event)
            while Pk ~= 0
                Pk = P(find(t==t_qd_s(i,j))-k);
                k = k+1;
            end
            t_p_s(i,j) = t(find(t==t_qd_s(i,j))-k+2);
            P_event = P(find(t==t_p_s(i,j)):find(t==t_p_e(i,j)));
            k = find(t==t_qd_s(i,j))-find(P(1:find(t==t_qd_s(i,j))-(k+1))~=0,1,'last');
            Pk = P(find(t==t_qd_s(i,j))-k);
        end
        
        % 5. calculate effective precipitation 
        % initial value for infiltration loss
        fi = 0;
        % infiltration loss
        sPloss = sum(P_event(fi >= P_event)) + sum(P_event > fi)*fi;
        % effective precipitation
        sPeff = sum(P_event) - sPloss;
        % sign and value of estimation error
        e = sPeff - R;
        s0 = sign(e); s = s0;
        % change increment for inf. loss
        dfi = 1;
        % initial array of loss and effective precipitation
        Ploss = zeros(size(P_event)); Peff = P_event;
        % iterative solver
        while abs(e) > R/10
            if s0 == s
                fi = fi+s*dfi;
            else
                dfi = dfi/10;
                fi = fi+s*dfi;
                s0 = s;
            end
            Ploss(fi >= P_event) = P_event(fi >= P_event);
            Ploss(P_event >= fi) = fi;
            Peff = P_event-Ploss;
            sPeff = sum(Peff);
            e = sPeff - R;
            s = sign(e);
        end
        FI(i,j) = fi;
        
        % 6. calculate the center of mass of effective precipitation
        t_p_event = hours(t(find(t==t_p_s(i,j)):find(t==t_p_e(i,j)))-t_p_s(i,j));
        P_com(i,j) = sum(Peff.*(Peff./2))/sum(Peff)+max(Ploss);
        t_P_com(i,j) = t_p_s(i,j) + duration([sum(Peff.*t_p_event)/sum(Peff) 0 0]);
        
        % 7. get end of effective precipitation
        t_Peff_e = t_p_s(i,j) + duration([find(Peff~=0,1,'last') 0 0]);
        
        % 8. calculate the center of mass of direct runoff
        t_qd_event = hours(t(find(t==t_qd_s(i,j)):find(t==t_qd_e(i,j)))-t_qd_s(i,j));
        qd_event = Qd(find(t==t_qd_s(i,j)):find(t==t_qd_e(i,j)));
        qd_com(i,j) = sum(qd_event.*(Qb(find(t==t_qd_s(i,j)):find(t==t_qd_e(i,j))) + qd_event./2))/sum(qd_event);
        t_qd_com(i,j) = t_qd_s(i,j) + duration([sum(qd_event.*t_qd_event)/sum(qd_event) 0 0]);
        
        % 9. Tc, Te, Tl
        Tc(i,j) = hours(t_qd_e(i,j) - t_Peff_e);
        Te(i,j) = hours(tpeak(j) - t_p_s(i,j));
        Tl(i,j) = hours(t_qd_com(i,j) - t_P_com(i,j));
        
        % 10. r = Tc/Te
        r(i,j) = Tc(i,j)/Te(i,j);
    end
    
end

% 11. final results
r_med = median(r,2,'omitnan');
beta_0 = beta(find(min(abs(r_med-2.25))==abs(r_med-2.25),1));
Qb_0 = qsep(Q,beta_0);
Tc_0 = Tc(beta==beta_0,:);
Te_0 = Te(beta==beta_0,:);
Tl_0 = Tl(beta==beta_0,:);

% 12. plotting
f = figure;
if plot_mode == 1
    for j = 1:numel(tpeak)
        id_peak = find(abs(hours(t-tpeak(j))) == min(abs(hours(t-tpeak(j)))),1,'first');
        
        clear Ploss
        fi = FI(beta==beta_0,j);
        P_event = P(find(t==t_p_s(beta==beta_0,j)):find(t==t_p_e(beta==beta_0,j)));
        Ploss(fi >= P_event) = P_event(fi >= P_event);
        Ploss(P_event >= fi) = fi;

        subplot(2,ceil(numel(tpeak)/2),j)
        plot(t,Qb_0,'k--','color',[.5 .5 .5],'linewidth',1.5)
        hold on
        plot(t,Q,'k.-','color',[.5 .5 .5])
        plot(t(find(t==t_qd_s(beta==beta_0,j)):find(t==t_qd_e(beta==beta_0,j))),...
            Q(find(t==t_qd_s(beta==beta_0,j)):find(t==t_qd_e(beta==beta_0,j))),'r.-')
        plot(t(id_peak),Q(id_peak),'ksq','Linewidth',1.5)
        plot([t_qd_s(beta==beta_0,j) t_qd_e(beta==beta_0,j)],...
            [Q(t==t_qd_s(beta==beta_0,j)) Q(t==t_qd_e(beta==beta_0,j))],'kx','linewidth',1.5)
        plot(t_qd_com(beta==beta_0,j),qd_com(beta==beta_0,j),'ko','Linewidth',1.5)
        ylim([0 2*Q(id_peak)])
        ylabel('Q [m^3/s]')

        yyaxis right
        bar(t,P,1,'FaceColor',[.5 .5 1],'EdgeColor','none')
        bar(t(find(t==t_p_s(beta==beta_0,j)):find(t==t_p_e(beta==beta_0,j))),...
            P(find(t==t_p_s(beta==beta_0,j)):find(t==t_p_e(beta==beta_0,j))),1,'r','EdgeColor','none')
        bar(t(find(t==t_p_s(beta==beta_0,j)):find(t==t_p_e(beta==beta_0,j))),Ploss,...
            1,'EdgeColor','none','FaceColor',[.5 .5 .5])
        plot([t_p_s(beta==beta_0,j) t_p_e(beta==beta_0,j)],[0 0],'kx','linewidth',1.5)
        plot(t_P_com(beta==beta_0,j),P_com(beta==beta_0,j),'ko','Linewidth',1.5)
        set(gca,'YDir','reverse')
        set(gca,'YColor','k')
        ylim([0 2*max(P(find(t==min([t_qd_s(beta==beta_0,j), t_p_s(beta==beta_0,j)])-duration([24 0 0])):...
            find(t==t_qd_e(beta==beta_0,j))))])
        ylabel('P [mm]')

        title(['T_c = ',num2str(round(Tc_0(j),1)),'; T_e = ',num2str(round(Te_0(j),1)),...
            '; T_L = ',num2str(round(Tl_0(j),1)),' [hr]'])
        xlim([min([t_qd_s(beta==beta_0,j), t_p_s(beta==beta_0,j)])-duration([24 0 0])...
            t_qd_e(beta==beta_0,j)+duration([24 0 0])])
        clear P_event Ploss
    end
    set(gcf,'position',[0,200,350*ceil(numel(tpeak)/2),600])
end

end

