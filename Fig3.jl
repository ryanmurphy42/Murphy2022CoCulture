# Murphy et al. (2022) - CoCulture Spheroids - Plots for Figure 3 (spheroid growth WM983B)

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
isdir(pwd() * "\\Fig3\\") || mkdir(pwd() * "\\Fig3") # make folder to save figures if doesnt already exist
filepath_save = [pwd() * "\\Fig3\\"] # location to save figures "\\Fig2\\"] # location to save figures

#######################################################################################
## Load data
data_combined_v2 = CSV.read(pwd() * "\\Fig3data.csv", DataFrame);

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

    unique_DaysSinceSeeding = (0:6:240)./24

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
        global plot1 = plot(unique_DaysSinceSeeding, mean_DaysSinceSeeding,yerr= std_DaysSinceSeeding,xlim=(3,10), ylim=(0,350), xlab=L"t \mathrm{ \ [days]}",ylab=L"R(t) \ [\mu \mathrm{m}]",xticks = [3,4,5,6,7,8,9,10],yticks=[0,100,200,300],titlefont=fnt, guidefont=fnt, tickfont=fnt,markersize = 5, lw=3,label=plot_Conditions[j]) #markerstrokecolor=:auto)
    else
        plot1 = plot!(unique_DaysSinceSeeding, mean_DaysSinceSeeding,yerr= std_DaysSinceSeeding,xlim=(3,10), ylim=(0,350), xlab=L"t \mathrm{ \ [days]}",ylab=L"R(t) \ [\mu \mathrm{m}]",xticks = [3,4,5,6,7,8,9,10],yticks=[0,100,200,300],titlefont=fnt, guidefont=fnt, tickfont=fnt,markersize = 5, lw=3,label=plot_Conditions[j],legend=:bottomright) #,markerstrokecolor=:auto
    end
end
display(plot1)
savefig(plot1,filepath_save[1] * "Fig1" * ".pdf")



#######################################################################################
### Loop through each condition for MLE, profiles, and bounds of confidence interval
## Parameter estimation using a linear model R(t) = R(3) + m*t

plot_Conditions =["M100:F0", "M75:F25", "M50:F50","M25:F75","M0:F100"]
for j=1:5

    # filter the data for condition j
    data_plot_condition = copy(data_combined_v2);
    function filter_plot_condition(Condition)::Bool
        interesting = Condition == plot_Conditions[j]
    end
    filter!([:Condition] => filter_plot_condition, data_plot_condition)

    global t1 = data_plot_condition[:,:DaysSinceSeeding] -3*ones(length(data_plot_condition[:,:DaysSinceSeeding]),1) # measurement times
    global data = data_plot_condition[:,:Radius] # radial measurements


    a=zeros(3)
    # define first guesses for MLE
    RR3=160;
    MM=0.1;
    SDD = 1.0;
    TH=-1.921; #95% confidence interval threshold

    function linear(t1,a)
        # linear model
        ylinear = a[1].*ones(length(t1),1) + a[1].*a[2].*t1;
        return(ylinear)
    end

    function error(data,a)
        # error model
        y=zeros(length(t1))
        y=linear(t1,a)
        e=0;
        dist=Normal(0,a[3]);
        e=loglikelihood(dist,data-y)
        ee=sum(e)
        return ee
    end

    function fun(a)
        return error(data,a)
    end

    function optimise(fun,θ₀,lb,ub)
        # optimisation function
        tomax = (θ,∂θ) -> fun(θ)
        opt = Opt(:LN_NELDERMEAD,length(θ₀))
        opt.max_objective = tomax
        opt.lower_bounds = lb       # Lower bound
        opt.upper_bounds = ub       # Upper bound
        opt.maxtime = 5
        res = optimize(opt,θ₀)
        return res[[2,1]]
    end

    #######################################################################################
    # MLE

    θG = [RR3,MM,SDD] # first guess for parameter estimation
    lb=[120.0,0.0,0.0]; # lower bound for parameter estimation
    ub=[200.0,0.2,20.0]; # upper bound for parameter estimation
    
    # MLE optimsation
    (xopt,fopt)  = optimise(fun,θG,lb,ub)
   
    # storing MLE 
    global fmle=fopt
    global R3mle=xopt[1]
    global Mmle=xopt[2]
    global Sdmle=xopt[3]

    # plot model simulated at MLE and save
    ymle = linear(t1,xopt);
    t1_smooth = LinRange(0,maximum(t1),10001)
    ymle_smooth = linear(t1_smooth,xopt);
    p1=scatter(t1 + 3*ones(length(t1),1),data,markersize = 3,markercolor=:black)
    p1=plot!(t1_smooth+3*ones(length(t1_smooth),1),ymle_smooth,xlab=L"t \mathrm{ \ [days]}",ylab=L"R(t) \ [\mu \mathrm{m}]",legend=false,xlims=(3,10),ylims=(0,350),xticks=[3,4,5,6,7,8,9,10],yticks=[0,100,200,300],lw=4,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=:orange)
    savefig(p1,filepath_save[1] * "Fig1" * replace(plot_Conditions[j],":" => "" )   * ".pdf")

    
    # #######################################################################################
    # # Profiling

    nptss=20

    #Profile MM
    Mmin=0.07
    if j==5
        Mmin=0.001
    end
    Mmax=0.1
    Mrange_lower=reverse(LinRange(Mmin,Mmle,nptss))
    Mrange_upper=LinRange(Mmle + (Mmax-Mmle)/nptss,Mmax,nptss)

    nrange_lower=zeros(2,nptss)
    llM_lower=zeros(nptss)
    nllM_lower=zeros(nptss)
    predict_M_lower=zeros(length(t1_smooth),nptss)

    nrange_upper=zeros(2,nptss)
    llM_upper=zeros(nptss)
    nllM_upper=zeros(nptss)
    predict_M_upper=zeros(length(t1_smooth),nptss)

    # start at mle and increase parameter (upper)
    for i in 1:nptss
        function fun1(aa)
            return error(data,[aa[1],Mrange_upper[i],aa[2]])
        end

        lb1=[lb[1],lb[3]];
        ub1=[ub[1],ub[3]];

        if i==1
            local θG1=[R3mle,Sdmle]
        elseif i==2
            # zero order approximation
            local θG1=nrange_upper[:,i-1]
        elseif i > 2
            # first order approximation
            local θG1= nrange_upper[:,i-1] + ((Mrange_upper[i]-Mrange_upper[i-1])./(Mrange_upper[i-1]-Mrange_upper[i-2]))*(nrange_upper[:,i-1]-nrange_upper[:,i-2])
            # if first order approximation is outside of lower bounds or upper bounds use zero order approximation
            if (sum(θG1.< lb1) + sum(θG1 .> ub1)) > 0
                local θG1=nrange_upper[:,i-1]
            end

        end
        error(data,[θG1[1],Mrange_upper[i],θG1[2]])


        local (xo,fo)=optimise(fun1,θG1,lb1,ub1)
        nrange_upper[:,i]=xo[:]
        llM_upper[i]=fo[1]
        if fo > fmle
            global fmle = fo
            global R3le=nrange_upper[1,i]
            global Mmle=Mrange_upper[i]
            global Sdmle=nrange_upper[2,i]
        end
    end

    # start at mle and decrease parameter (lower)
    for i in 1:nptss
        function fun1a(aa)
            return error(data,[aa[1],Mrange_lower[i],aa[2]])
        end

        lb1=[lb[1],lb[3]];
        ub1=[ub[1],ub[3]];

        if i==1
            local θG1=[R3mle,Sdmle]
        elseif i==2
            # zero order approximation
            local θG1=nrange_lower[:,i-1]
        elseif i > 2
            # first order approximation
            local θG1= nrange_lower[:,i-1] + ((Mrange_lower[i]-Mrange_lower[i-1])./(Mrange_lower[i-1]-Mrange_lower[i-2]))*(nrange_lower[:,i-1]-nrange_lower[:,i-2])
            # if first order approximation is outside of lower bounds or upper bounds use zero order approximation
            if (sum(θG1.< lb1) + sum(θG1 .> ub1)) > 0
                local θG1=nrange_lower[:,i-1]
            end

        end

        local (xo,fo)=optimise(fun1a,θG1,lb1,ub1)
        nrange_lower[:,i]=xo[:]
        llM_lower[i]=fo[1]
        if fo > fmle
            global fmle = fo
            global R3le=nrange_lower[1,i]
            global Mmle=Mrange_lower[i]
            global Sdmle=nrange_lower[2,i]
        end
    end

    # combine the lower and upper
    Mrange = [reverse(Mrange_lower); Mrange_upper]
    nrange = [reverse(nrange_lower); nrange_upper ]
    llM = [reverse(llM_lower); llM_upper]

    nllM=llM.-maximum(llM);



    ## Profile R3
    R3min=140
    R3max=200
    R3range_lower=reverse(LinRange(R3min,R3mle,nptss))
    R3range_upper=LinRange(R3mle + (R3max-R3mle)/nptss,R3max,nptss)

    nrange_lower=zeros(2,nptss)
    llR3_lower=zeros(nptss)
    nllR3_lower=zeros(nptss)
    predict_R3_lower=zeros(length(t1_smooth),nptss)

    nrange_upper=zeros(2,nptss)
    llR3_upper=zeros(nptss)
    nllR3_upper=zeros(nptss)
    predict_R3_upper=zeros(length(t1_smooth),nptss)

    cur_colors = palette(:default)


    # start at mle and increase parameter (upper)
    for i in 1:nptss
        function fun1(aa)
            return error(data,[R3range_upper[i],aa[1],aa[2]])
        end

        lb1=[lb[2],lb[3]];
        ub1=[ub[3],ub[3]];

        if i==1
            local θG1=[Mmle,Sdmle]
        elseif i==2
            # zero order approximation
            local θG1=nrange_upper[:,i-1]
        elseif i > 2
            # first order approximation
            local θG1= nrange_upper[:,i-1] + ((R3range_upper[i]-R3range_upper[i-1])./(R3range_upper[i-1]-R3range_upper[i-2]))*(nrange_upper[:,i-1]-nrange_upper[:,i-2])
            # if first order approximation is outside of lower bounds or upper bounds use zero order approximation
            if (sum(θG1.< lb1) + sum(θG1 .> ub1)) > 0
                local θG1=nrange_upper[:,i-1]
            end

        end

        local (xo,fo)=optimise(fun1,θG1,lb1,ub1)
        nrange_upper[:,i]=xo[:]
        llR3_upper[i]=fo[1]
        
        if fo > fmle
            global fmle = fo
            global R3le=R3range_upper[i]
            global Mmle=nrange_upper[1,i]
            global Sdmle=nrange_upper[2,i]
        end
    end

    # start at mle and decrease parameter (lower)
    for i in 1:nptss
        function fun1a(aa)
            return error(data,[R3range_lower[i],aa[1],aa[2]])
        end

        lb1=[lb[2],lb[3]];
        ub1=[ub[2],ub[3]];

        if i==1
            local θG1=[Mmle,Sdmle]
        elseif i==2
            # zero order approximation
            local θG1=nrange_lower[:,i-1]
        elseif i > 2
            # first order approximation
            local θG1= nrange_lower[:,i-1] + ((R3range_lower[i]-R3range_lower[i-1])./(R3range_lower[i-1]-R3range_lower[i-2]))*(nrange_lower[:,i-1]-nrange_lower[:,i-2])
            # if first order approximation is outside of lower bounds or upper bounds use zero order approximation
            if (sum(θG1.< lb1) + sum(θG1 .> ub1)) > 0
                local θG1=nrange_lower[:,i-1]
            end

        end

        local (xo,fo)=optimise(fun1a,θG1,lb1,ub1)
        nrange_lower[:,i]=xo[:]
        llR3_lower[i]=fo[1]

        if fo > fmle
            global fmle = fo
            global R3le=R3range_lower[i]
            global Mmle=nrange_lower[1,i]
            global Sdmle=nrange_lower[2,i]
        end
    end

    # combine the lower and upper
    R3range = [reverse(R3range_lower); R3range_upper]
    nrange = [reverse(nrange_lower); nrange_upper ]
    llR3 = [reverse(llR3_lower); llR3_upper]

    nllR3=llR3.-maximum(llR3);


     # interpolate for smoother profile likelihoods
     interp_nptss= 1001;

     interp_points_Mrange =  LinRange(Mmin,Mmax,interp_nptss)
     interp_M = LinearInterpolation(Mrange,nllM)
     interp_nllM = interp_M(interp_points_Mrange)

       s1=plot(interp_points_Mrange,interp_nllM,xlim=(0.06,0.1),ylim=(-4,0.1),xticks=[0.05, 0.06,0.07,0.08,0.09, 0.1, 0.2],yticks=[-3,-2,-1,0],xlab=L"\lambda \ [\mathrm{day}^{-1}]",ylab=L"\hat{\ell}_{p}",legend=false,lw=5,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=:deepskyblue3)
       s1=hline!([-1.92],lw=2,linecolor=:black,linestyle=:dot)
       s1=vline!([Mmle],lw=3,linecolor=:red)

        display(s1)
        savefig(s1,filepath_save[1] * "Figs1" * replace(plot_Conditions[j],":" => "" )   * ".pdf")


        # plot profile likelihoods on same figures
        if j==1
        global s1_combined = hline([-1.92],lw=2,linecolor=:black,linestyle=:dot)
        end
        s1_combined=plot!(s1_combined,interp_points_Mrange,interp_nllM,xlim=(0.07,0.1),ylim=(-4,0.1),xticks=[0.05, 0.06,0.07,0.08,0.09, 0.1, 0.2],yticks=[-3,-2,-1,0],xlab=L"\lambda \ [\mathrm{day}^{-1}]",ylab=L"\hat{\ell}_{p}",legend=false,lw=2,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=cur_colors[j])
        s1_combined=vline!(s1_combined,[Mmle],lw=2,linecolor=cur_colors[j],linestyle=:dash)

        display(s1_combined)
        savefig(s1_combined,filepath_save[1] * "Figs1_combined" * replace(plot_Conditions[j],":" => "" )   * ".pdf")



    interp_points_R3range =  LinRange(R3min,R3max,interp_nptss)
    interp_R3 = LinearInterpolation(R3range,nllR3)
    interp_nllR3 = interp_R3(interp_points_R3range)

    s2=plot(interp_points_R3range,interp_nllR3,xlim=(100,200),ylim=(-4,0.1),xticks=[100,150,200],yticks=[-3,-2,-1,0],xlab=L"R(3) \ [\mu\mathrm{m}]",ylab=L"\hat{\ell}_{p}",legend=false,lw=5,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=:deepskyblue3)
    s2=hline!([-1.92],lw=2,linecolor=:black,linestyle=:dot)

     display(s2)
     savefig(s2,filepath_save[1] * "Figs2" * replace(plot_Conditions[j],":" => "" )   * ".pdf")

    if j==1
        global s2_combined=hline([-1.92],lw=2,linecolor=:black,linestyle=:dot)
    end
     s2_combined=plot!(s2_combined,interp_points_R3range,interp_nllR3,xlim=(150,190),ylim=(-4,0.1),xticks=[150,160,170,180,190,200],yticks=[-3,-2,-1,0],xlab=L"R(3) \ [\mu \mathrm{m}]",ylab=L"\hat{\ell}_{p}",legend=false,lw=2,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=cur_colors[j])
     s2_combined=vline!(s2_combined,[R3mle],lw=2,linecolor=cur_colors[j],linestyle=:dash)

     display(s2_combined)
     savefig(s2_combined,filepath_save[1] * "Figs2_combined" * replace(plot_Conditions[j],":" => "" )   * ".pdf")

       function fun_interpCI(mle,interp_points_range,interp_nll,TH)
           # find bounds of CI

           range_minus_mle = interp_points_range - mle*ones(length(interp_points_range),1)
           abs_range_minus_mle = broadcast(abs, range_minus_mle)
           findmin_mle = findmin(abs_range_minus_mle)
           # find closet value to CI threshold intercept
           value_minus_threshold = interp_nll - TH*ones(length(interp_nll),1)
           abs_value_minus_threshold = broadcast(abs, value_minus_threshold)
           lb_CI_tmp = findmin(abs_value_minus_threshold[1:findmin_mle[2][1]])
           ub_CI_tmp = findmin(abs_value_minus_threshold[findmin_mle[2][1]:length(abs_value_minus_threshold)])
           lb_CI = interp_points_range[lb_CI_tmp[2][1]]
           ub_CI = interp_points_range[findmin_mle[2][1]-1 + ub_CI_tmp[2][1]]

           return lb_CI,ub_CI
       end

       # M
       (lb_CI_M,ub_CI_M) = fun_interpCI(Mmle,interp_points_Mrange,interp_nllM,TH)
       println(round(lb_CI_M; digits = 2))
       println(round(ub_CI_M; digits = 2))

        # R3
        (lb_CI_R3,ub_CI_R3) = fun_interpCI(R3mle,interp_points_R3range,interp_nllR3,TH)
        println(round(lb_CI_R3; digits = 2))
        println(round(ub_CI_R3; digits = 2))


        # Export MLE and bounds to csv (one file for all data) -- -MLE ONLY

        println(@isdefined(df_MLEBoundsAll) == 0)

        if @isdefined(df_MLEBoundsAll) == 0
            println("not defined")
            global df_MLEBoundsAll = DataFrame(Condition = replace(plot_Conditions[j],":" => "" ), lambdamle=Mmle,lb_CI_lambda=lb_CI_lambda,ub_CI_M=ub_CI_M,R3mle=R3mle,lb_CI_R3=lb_CI_R3,ub_CI_R3=ub_CI_R3)
        else
            println("DEFINED")
            global df_MLEBoundsAll_thisrow = DataFrame(Condition = replace(plot_Conditions[j],":" => "" ), lambdamle=Mmle,lb_CI_lambda=lb_CI_M,ub_CI_lambda=ub_CI_M,R3mle=R3mle,lb_CI_R3=lb_CI_R3,ub_CI_R3=ub_CI_R3)
            append!(df_MLEBoundsAll,df_MLEBoundsAll_thisrow)
        end

        CSV.write(filepath_save[1] * "MLEBoundsALL.csv", df_MLEBoundsAll)


end
