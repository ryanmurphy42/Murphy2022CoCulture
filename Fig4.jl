# Murphy et al. (2022) - CoCulture Spheroids - Plots for Figure 4 (biphasic model 1205LU)

#######################################################################################
## Initialisation including packages, plot options, filepaths to save outputs

using Plots
using LinearAlgebra
using NLopt
using .Threads 
using Interpolations
using Distributions
using Roots
using LaTeXStrings
using CSV
using DataFrames

pyplot() # plot options
fnt = Plots.font("sans-serif", 20) # plot options
global cur_colors = palette(:default) # plot options
isdir(pwd() * "\\Fig4\\") || mkdir(pwd() * "\\Fig4") # make folder to save figures if doesnt already exist
filepath_save = [pwd() * "\\Fig4\\"] # location to save figures "\\Fig4\\"] # location to save figures

#######################################################################################
## Load data
data_combined_v2 = CSV.read(pwd() * "\\Fig4data.csv", DataFrame);

#######################################################################################
### Plot data - mean plusminus standard deviation

plot_Conditions =["M100:F0", "M75:F25", "M50:F50","M25:F75","M0:F100"]
for j=1:5
    # filter the data for condition j
    data_plot_condition = copy(data_combined_v2);
    function filter_plot_condition(Condition)::Bool
        interesting = Condition == plot_Conditions[j]
    end
    filter!([:Condition] => filter_plot_condition, data_plot_condition)

    unique_DaysSinceSeeding = unique(data_combined_v2[:,:DaysSinceSeeding])

    # calculate mean and standard deviation of data at each timepoint
    mean_DaysSinceSeeding = zeros(length(unique_DaysSinceSeeding),1);
    std_DaysSinceSeeding = zeros(length(unique_DaysSinceSeeding),1);
    for i1 = 1:length(unique_DaysSinceSeeding)
        data_plot_condition_pertimepoint = copy(data_plot_condition);
        function filter_plot_condition_pertimepoint(DaysSinceSeeding)::Bool
            interesting = DaysSinceSeeding == unique_DaysSinceSeeding[i1]
        end
        filter!([:DaysSinceSeeding] => filter_plot_condition_pertimepoint, data_plot_condition_pertimepoint)
        mean_DaysSinceSeeding[i1] = mean(data_plot_condition_pertimepoint[:,:Radius])
        std_DaysSinceSeeding[i1] = std(data_plot_condition_pertimepoint[:,:Radius])
    end
    # plot data
    if j==1
        global plot1 = plot(unique_DaysSinceSeeding, mean_DaysSinceSeeding,yerr= std_DaysSinceSeeding,xlim=(0,10), ylim=(0,700), xlab=L"t \mathrm{ \ [days]}",ylab=L"R(t) \ [\mu \mathrm{m}]",xticks = [0,2,4,6,8,10],yticks=[0,100,200,300,400,500,600],titlefont=fnt, guidefont=fnt, tickfont=fnt,markersize = 5, lw=3,label=plot_Conditions[j]) #markerstrokecolor=:auto)
        display(plot1)
    else
        global plot1 = plot!(plot1,unique_DaysSinceSeeding, mean_DaysSinceSeeding,yerr= std_DaysSinceSeeding,xlim=(0,10), ylim=(0,700), xlab=L"t \mathrm{ \ [days]}",ylab=L"R(t) \ [\mu \mathrm{m}]",xticks = [0,2,4,6,8,10],yticks=[0,100,200,300,400,500,600],titlefont=fnt, guidefont=fnt, tickfont=fnt,markersize = 5, lw=3,label=plot_Conditions[j],legend=:bottomright) #,markerstrokecolor=:auto
        display(plot1)
    end
end

savefig(plot1,filepath_save[1] * "Fig1" * ".pdf")

#######################################################################################
### Loop through each condition for MLE, profiles, and bounds of confidence interval


plot_Conditions =["M100:F0", "M75:F25", "M50:F50","M25:F75","M0:F100"]
for j=1:4

        # filter the data for condition j
        data_plot_condition = copy(data_combined_v2);
        function filter_plot_condition(Condition)::Bool
            interesting = Condition == plot_Conditions[j]
        end
        filter!([:Condition] => filter_plot_condition, data_plot_condition)

        
        t1 = data_plot_condition[:,:DaysSinceSeeding] -(day_shift)*ones(length(data_plot_condition[:,:DaysSinceSeeding]),1); # measurement times
        data = data_plot_condition[:,:Radius]; # radial measurements

        # set mean of first measurements as initial guess for initial size
        data_first_measurement = copy(data_plot_condition);
        function filter_first_measurement(DaysSinceSeeding)::Bool
            interesting = DaysSinceSeeding == day_shift
        end
        filter!([:DaysSinceSeeding] => filter_first_measurement, data_first_measurement)

        # set first guess for KK - maximum size at final time
        data_last_measurement = copy(data_plot_condition);
        function filter_last_measurement(DaysSinceSeeding)::Bool
            interesting = DaysSinceSeeding == maximum(data_last_measurement[:,:DaysSinceSeeding])
        end
        filter!([:DaysSinceSeeding] => filter_last_measurement, data_last_measurement)

    a=zeros(7)
    # initial guesses for MLE search
    rr=0.1;
    dd=1;
    KK=round(maximum(data_last_measurement[:,:Radius]),digits=2);
    JJ=round(minimum(data_plot_condition[:,:Radius]),digits=2);
    CC0=round(mean(data_first_measurement[:,:Radius]),digits=2);
    TT=2.0;
    Sd=10.0
    TH=-1.921; #95% confidence interval threshold

    function logistic_delayed(t1,a)
        # function to compute the biphasic logistic model
        r=a[1];
        d=a[2]
        K=a[3];
        J=a[4];
        C0=a[5];
        TT=a[6];
        dt=maximum(t1)/10000;
        t1mesh=0:dt:TT;
        t2mesh=TT+dt:dt:maximum(t1)+5*dt;
        tmesh=vcat(t1mesh,t2mesh);
        C=zeros(length(tmesh));
        f1 = (c,d) -> d*c*(1-c/J);
        f2 = (c,r,K) -> r*c*(1-c/K);
        C[1]=C0;
            for i in 2:length(t1mesh)
                ctemp = C[i-1]+ f1(C[i-1],d)*dt;
                C[i] =C[i-1]+ 0.5*(f1(ctemp,d)+f1(C[i-1],d))*dt;
            end
            for i in length(t1mesh)+1:length(tmesh)
                ctemp = C[i-1]+ f2(C[i-1],r,K)*dt;
                C[i] =C[i-1]+ 0.5*(f2(ctemp,r,K)+f2(C[i-1],r,K))*dt;
            end 
            f=LinearInterpolation(tmesh,C)
            interp=zeros(length(t1))
            interp=f(t1[1:length(t1)])
        return interp
    end


    function error(data,a)
        # error model
        y=zeros(length(t1))
        y=logistic_delayed(t1,a);
        e=0;
        dist=Normal(0,a[7]);
        e=loglikelihood(dist,data-y) 
        ee=sum(e)
        return ee
    end


    function fun(a)
        return error(data,a)
    end


    function optimise(fun,θ₀,lb,ub;
        # optimisation function
        dv = false,
        method = dv ? :LD_LBFGS : :LN_BOBYQA)

        if dv || String(method)[2] == 'D'
        tomax = fun
        else
            tomax = (θ,∂θ) -> fun(θ)
        end

        opt = Opt(method,length(θ₀))
        opt.max_objective = tomax
        opt.lower_bounds = lb       # Lower bound
        opt.upper_bounds = ub       # Upper bound
        opt.local_optimizer = Opt(:LN_NELDERMEAD, length(θ₀))
        opt.maxtime = 30.0; # maximum time in seconds
        res = optimize(opt,θ₀)
        return res[[2,1]]
    end


    #######################################################################################
    # MLE

    
    θG = [rr,dd,KK,JJ,CC0,TT,Sd] # first guess
    # lower and upper bounds for parameter estimation
    r_lb = 0.00001
    r_ub = 2.0
    d_lb = 0.00001
    d_ub = 2.0
    K_lb = 100.0
    K_ub = 600.0
    J_lb = 100
    J_ub = 400
    C0_lb = 100.0
    C0_ub = 700.0
    T_lb = 0.0
    T_ub = 5.0
    V_lb = 0.0
    V_ub = 20.0

    lb=[r_lb,d_lb,K_lb,J_lb,C0_lb,T_lb,V_lb];
    ub=[r_ub,d_ub,K_ub,J_ub,C0_ub,T_ub,V_ub];

    # MLE optimisation
    (xopt,fopt)  = optimise(fun,θG,lb,ub)

    # storing MLE 
    global fmle=fopt
    global rmle=xopt[1]
    global dmle=xopt[2]
    global Kmle=xopt[3]
    global Jmle=xopt[4]
    global C0mle=xopt[5]
    global Tmle=xopt[6]
    global Sdmle=xopt[7]

    # plot model simulated at MLE and save
    ymle = logistic_delayed(t1,xopt);
    t1_smooth = LinRange(0,maximum(t1),10001)
    ymle_smooth = logistic_delayed(t1_smooth,xopt);
    p1=scatter(t1+day_shift*ones(length(t1),1),data,xlims=(0,10),ylims=(0,700),legend=false,markersize = 3,markercolor=:black)
    p1=plot!(t1_smooth + day_shift*ones(length(t1_smooth),1),ymle_smooth,xticks=[0,2,4,6,8,10],yticks=[0,100,200,300,400,500,600],xlab=L"t",ylab=L"R(t)",lw=4,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=:orange)
    display(p1)
    savefig(p1,filepath_save[1] * "Figp1" * replace(plot_Conditions[j],":" => "" ) * ".pdf")


    # #######################################################################################
    # # Profiling
    nptss=40

    #Profile T
    Tmin=1-day_shift
    Tmax=3-day_shift
    Trange_lower=reverse(LinRange(Tmin,Tmle,nptss))
    Trange_upper=LinRange(Tmle + (Tmax-Tmle)/nptss,Tmax,nptss)

    nrange_lower=zeros(6,nptss)
    llT_lower=zeros(nptss)
    nllT_lower=zeros(nptss)
    predict_T_lower=zeros(length(t1_smooth),nptss)

    nrange_upper=zeros(6,nptss)
    llT_upper=zeros(nptss)
    nllT_upper=zeros(nptss)
    predict_T_upper=zeros(length(t1_smooth),nptss)

    # start at mle and increase parameter (upper)
    for i in 1:nptss
        function fun4(aa)
            return error(data,[aa[1],aa[2],aa[3],aa[4],aa[5],Trange_upper[i],aa[6]])
        end

        lb1=[r_lb,d_lb,K_lb,J_lb,C0_lb,V_lb];
        ub1=[r_ub,d_ub,K_ub,J_ub,C0_ub,V_ub];
            
        if i==1
            local θG1=[rmle,dmle,Kmle,Jmle,C0mle,Sdmle]
        elseif i==2
            # zero order approximation
            local θG1=nrange_upper[:,i-1]
            # if zero order approximation outside of bounds use MLE
            if (sum(θG1.< lb1) + sum(θG1 .> ub1)) > 0
                local θG1=[rmle,dmle,Kmle,Jmle,C0mle,Sdmle]
            end
        elseif i > 2
            # first order approximation
            local θG1= nrange_upper[:,i-1] + ((Trange_upper[i]-Trange_upper[i-1])./(Trange_upper[i-1]-Trange_upper[i-2]))*(nrange_upper[:,i-1]-nrange_upper[:,i-2])
            # if first order approximation is outside of lower bounds or upper bounds use zero order approximation
            if (sum(θG1.< lb1) + sum(θG1 .> ub1)) > 0
                #local θG1=nrange_upper[:,i-1]
                local θG1=[rmle,dmle,Kmle,Jmle,C0mle,Sdmle]
            end
            # if zero order approximation outside of bounds use MLE
            if (sum(θG1.< lb1) + sum(θG1 .> ub1)) > 0
                local θG1=[rmle,dmle,Kmle,Jmle,C0mle,Sdmle]
            end
        end

        local (xo,fo)=optimise(fun4,θG1,lb1,ub1)
        nrange_upper[:,i]=xo[:]
        llT_upper[i]=fo[1]
        predict_T_upper[:,i]=logistic_delayed(t1_smooth,[nrange_upper[1,i],nrange_upper[2,i],nrange_upper[3,i],nrange_upper[4,i],nrange_upper[5,i],Trange_upper[i],nrange_upper[6,i]])

        if fo > fmle
            global fmle = fo
            global rmle=nrange_upper[1,i]
            global dmle=nrange_upper[2,i]
            global Kmle=nrange_upper[3,i]
            global Jmle=nrange_upper[4,i]
            global C0mle=nrange_upper[5,i]
            global Tmle=Trange_upper[i]
            global Sdmle=nrange_upper[6,i]
        end
    end

    # start at mle and decrease parameter (lower)
    for i in 1:nptss
        function fun4a(aa)
            return error(data,[aa[1],aa[2],aa[3],aa[4],aa[5],Trange_lower[i],aa[6]])
        end

        lb1=[r_lb,d_lb,K_lb,J_lb,C0_lb,V_lb];
        ub1=[r_ub,d_ub,K_ub,J_ub,C0_ub,V_ub];
        
        if i==1
            local θG1=[rmle,dmle,Kmle,Jmle,C0mle,Sdmle]
        elseif i==2
            # zero order approximation
            local θG1=nrange_lower[:,i-1]
            # if zero order approximation outside of bounds use MLE
            if (sum(θG1.< lb1) + sum(θG1 .> ub1)) > 0
                local θG1=[rmle,dmle,Kmle,Jmle,C0mle,Sdmle]
            end
        elseif i > 2
            # first order approximation
            local θG1= nrange_lower[:,i-1] + ((Trange_lower[i]-Trange_lower[i-1])./(Trange_lower[i-1]-Trange_lower[i-2]))*(nrange_lower[:,i-1]-nrange_lower[:,i-2])
            # if first order approximation is outside of lower bounds or upper bounds use zero order approximation
            if (sum(θG1.< lb1) + sum(θG1 .> ub1)) > 0
                #local θG1=nrange_lower[:,i-1]
                local θG1=[rmle,dmle,Kmle,Jmle,C0mle,Sdmle]
            end
            # if zero order approximation outside of bounds use MLE
            if (sum(θG1.< lb1) + sum(θG1 .> ub1)) > 0
                local θG1=[rmle,dmle,Kmle,Jmle,C0mle,Sdmle]
            end
        end

        local (xo,fo)=optimise(fun4a,θG1,lb1,ub1)
        nrange_lower[:,i]=xo[:]
        llT_lower[i]=fo[1]
        predict_T_lower[:,i]=logistic_delayed(t1_smooth,[nrange_lower[1,i],nrange_lower[2,i],nrange_lower[3,i],nrange_lower[4,i],nrange_lower[5,i],Trange_lower[i],nrange_lower[6,i]])

        if fo > fmle
            global fmle = fo
            global rmle=nrange_lower[1,i]
            global dmle=nrange_lower[2,i]
            global Kmle=nrange_lower[3,i]
            global Jmle=nrange_lower[4,i]
            global C0mle=nrange_lower[5,i]
            global Tmle=Trange_lower[i]
            global Sdmle=nrange_lower[6,i]
        end
    end

    # combine the lower and upper
    Trange = [reverse(Trange_lower);Trange_upper]
    llT = [reverse(llT_lower); llT_upper] 
    predict_T = [reverse(predict_T_lower,dims=2) predict_T_upper]

    nllT=llT.-maximum(llT);

    upper_T=zeros(length(t1_smooth))
    lower_T=1000*ones(length(t1_smooth))

    for i in 1:(nptss*2)
        if nllT[i] >= TH
            for j in 1:length(t1_smooth)
                upper_T[j]=max(predict_T[j,i],upper_T[j])
                lower_T[j]=min(predict_T[j,i],lower_T[j])
            end
        end
    end

    # interpolate for smoother profile likelihoods
    interp_nptss= 1001;

    # T
    interp_points_Trange =  LinRange(Tmin,Tmax,interp_nptss)
    interp_T = LinearInterpolation(Trange,nllT)
    interp_nllT = interp_T(interp_points_Trange)

    s4=plot!(interp_points_Trange+day_shift*ones(length(interp_points_Trange),1),interp_nllT,xlim=(0,4),ylim=(-4,0.1),xticks=[0,1,2,3,4],yticks=[-3,-2,-1,0],xlab=L"T \mathrm{ \ [days]}",ylab=L"\hat{\ell}_{p}",legend=false,lw=5,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=:deepskyblue3)
    s4=hline!([-1.92],lw=2,linecolor=:black,linestyle=:dot)
    s4=vline!([Tmle+day_shift],lw=3,linecolor=:red)

    if j==1
        global s4_combined=plot(interp_points_Trange+day_shift*ones(length(interp_points_Trange),1),interp_nllT,xlim=(0,4),ylim=(-4,0.1),xticks=[0,1,2,3,4],yticks=[-3,-2,-1,0],xlab=L"T \mathrm{ \ [days]}",ylab=L"\hat{\ell}_{p}",legend=false,lw=5,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=cur_colors[j])
        global s4_combined=hline!([-1.92],lw=2,linecolor=:black,linestyle=:dot)
        global s4_combined=vline!([Tmle+day_shift],lw=3,linestyle=:dash,linecolor=cur_colors[j])
    else
        global s4_combined=plot!(s4_combined,interp_points_Trange+day_shift*ones(length(interp_points_Trange),1),interp_nllT,xlim=(0,4),ylim=(-4,0.1),xticks=[0,1,2,3,4],yticks=[-3,-2,-1,0],xlab=L"T \mathrm{ \ [days]}",ylab=L"\hat{\ell}_{p}",legend=false,lw=5,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=cur_colors[j])
        global s4_combined=hline!([-1.92],lw=2,linecolor=:black,linestyle=:dot)
        global s4_combined=vline!([Tmle+day_shift],lw=3,linestyle=:dash,linecolor=cur_colors[j])
    end

    # Recompute best fit
    t1_smooth_updated= LinRange(0,maximum(t1),10001)
    ymle_smooth_updated = logistic_delayed(t1_smooth_updated,[rmle;dmle;Kmle;Jmle;C0mle;Tmle;Sdmle]);
    p1_updated=scatter(t1+day_shift*ones(length(t1),1),data,xlims=(0,10),ylims=(0,700),legend=false,markersize = 3,markercolor=:black)
    p1_updated=plot!(t1_smooth_updated + day_shift*ones(length(t1_smooth_updated),1),ymle_smooth_updated,xticks=[0,2,4,6,8,10],yticks=[0,100,200,300,400,500,600],xlab=L"t",ylab=L"R(t)",lw=4,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=:orange)

    # Save figures
    display(p1_updated)
    savefig(p1_updated,filepath_save[1] * "Figp1_updated" * replace(plot_Conditions[j],":" => "" ) * ".pdf")

    display(s4)
    savefig(s4,filepath_save[1] * "Figs4" * replace(plot_Conditions[j],":" => "" ) * ".pdf")

    display(s4_combined)
    savefig(s4_combined,filepath_save[1] * "Figs4_combined" * replace(plot_Conditions[j],":" => "" ) * ".pdf")

    #######################################################################################
    # compute the bounds of confidence interval

    function fun_interpCI(mle,interp_points_range,interp_nll,TH)
        # find bounds of CI
        range_minus_mle = interp_points_range - mle*ones(length(interp_points_range),1)
        abs_range_minus_mle = broadcast(abs, range_minus_mle)
        findmin_mle = findmin(abs_range_minus_mle)

        # find closest value to CI threshold intercept
        value_minus_threshold = interp_nll - TH*ones(length(interp_nll),1)
        abs_value_minus_threshold = broadcast(abs, value_minus_threshold)
        lb_CI_tmp = findmin(abs_value_minus_threshold[1:findmin_mle[2][1]])
        ub_CI_tmp = findmin(abs_value_minus_threshold[findmin_mle[2][1]:length(abs_value_minus_threshold)])
        lb_CI = interp_points_range[lb_CI_tmp[2][1]]
        ub_CI = interp_points_range[findmin_mle[2][1]-1 + ub_CI_tmp[2][1]]

        return lb_CI,ub_CI
    end

    # T
    (lb_CI_T,ub_CI_T) = fun_interpCI(Tmle,interp_points_Trange,interp_nllT,TH)
    println(plot_Conditions[j])
    println(round(lb_CI_T + day_shift; digits = 2)) 
    println(round(ub_CI_T + day_shift; digits = 2))
    println(round(Tmle + day_shift; digits = 2))
    println("done " * string(j))

end
