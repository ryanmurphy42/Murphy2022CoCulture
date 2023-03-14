# Murphy et al. (2022) - CoCulture Spheroids - Plots for Figure 7 (reduced Greenspan spheroid growth WM983B)

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
using DifferentialEquations
using Parsers # for writeCSV to work

pyplot() # plot options
fnt = Plots.font("sans-serif", 20) # plot options
global cur_colors = palette(:default) # plot options
isdir(pwd() * "\\Fig7\\") || mkdir(pwd() * "\\Fig7") # make folder to save figures if doesnt already exist
filepath_save = [pwd() * "\\Fig7\\"] # location to save figures "\\Fig2\\"] # location to save figures


#######################################################################################
## Load data
data_all_v1 = CSV.read(pwd() * "\\Fig7data.csv", DataFrame);


#######################################################################################
### Plot data - mean plusminus standard deviation
 
plot_Conditions =["M100:F0", "M75:F25", "M50:F50","M25:F75"]
for j=1:4
     # filter the data for condition j
    data_plot_condition = copy(data_all_v1);
    function filter_plot_condition(Condition)::Bool
        interesting = Condition == plot_Conditions[j]
    end
    filter!([:Condition] => filter_plot_condition, data_plot_condition)

    unique_DaysSinceSeeding = [2,3,6,8,10];

     # calculate mean and standard deviation of data at each timepoint
    mean_DaysSinceSeeding = zeros(length(unique_DaysSinceSeeding),1);
    std_DaysSinceSeeding = zeros(length(unique_DaysSinceSeeding),1);

    meanNecrotic_DaysSinceSeeding = zeros(length(unique_DaysSinceSeeding),1);
    stdNecrotic_DaysSinceSeeding = zeros(length(unique_DaysSinceSeeding),1);
    for i1 = 1:length(unique_DaysSinceSeeding)
        data_plot_condition_pertimepoint = copy(data_plot_condition);
        function filter_plot_condition_pertimepoint(DaysSinceSeeding)::Bool
            interesting = DaysSinceSeeding == unique_DaysSinceSeeding[i1]
        end
        filter!([:DaysSinceSeeding] => filter_plot_condition_pertimepoint, data_plot_condition_pertimepoint)
        mean_DaysSinceSeeding[i1] = mean(data_plot_condition_pertimepoint[:,:Radius])
        std_DaysSinceSeeding[i1] = std(data_plot_condition_pertimepoint[:,:Radius])
        meanNecrotic_DaysSinceSeeding[i1] = mean(data_plot_condition_pertimepoint[:,:NecroticRadius])
        stdNecrotic_DaysSinceSeeding[i1] = std(data_plot_condition_pertimepoint[:,:NecroticRadius])
    end
     # plot data and save for each condition
    global plot1 = plot(unique_DaysSinceSeeding, mean_DaysSinceSeeding,yerr= std_DaysSinceSeeding,xlim=(0,10), ylim=(0,350), xlab=L"t \mathrm{ \ [days]}",ylab=L"R(t) \ [\mu \mathrm{m}]",xticks = [0,2,4,6,8,10],yticks=[0,100,200,300],titlefont=fnt, guidefont=fnt, tickfont=fnt,markersize = 5, lw=3,label=plot_Conditions[j]) #markerstrokecolor=:auto)
    plot1 = plot!(unique_DaysSinceSeeding, meanNecrotic_DaysSinceSeeding,yerr= stdNecrotic_DaysSinceSeeding,xlim=(0,10), ylim=(0,350), xlab=L"t \mathrm{ \ [days]}",ylab=L"R(t) \ [\mu \mathrm{m}]",xticks = [0,2,4,6,8,10],yticks=[0,100,200,300],titlefont=fnt, guidefont=fnt, tickfont=fnt,markersize = 5, lw=3,label=plot_Conditions[j]) #markerstrokecolor=:auto)
    display(plot1)
    savefig(plot1,filepath_save[1] * "Fig1" * "_" * replace(plot_Conditions[j],":" => "" ) * ".pdf")
end


#######################################################################################
### Plot data - scatter Plots

plot_Conditions =["M100:F0", "M75:F25", "M50:F50","M25:F75"]
for j=1:4
    # filter data for condition j
    data_plot_condition = copy(data_all_v1);
    function filter_plot_condition(Condition)::Bool
        interesting = Condition == plot_Conditions[j]
    end
    filter!([:Condition] => filter_plot_condition, data_plot_condition)
    # plot data and save for each condition
    plot1a = scatter(data_plot_condition[:,:DaysSinceSeeding], data_plot_condition[:,:Radius],xlim=(0,11), ylim=(0,350), xlab=L"t \mathrm{ \ [days]}",ylab=L"R(t) \ [\mu \mathrm{m}]",xticks = [0,2,4,6,8,10],yticks=[0,100,200,300],titlefont=fnt, guidefont=fnt, tickfont=fnt,markersize = 5, lw=3,legend = false,color=:green2,mode="markers",markerstrokecolor=:green2)
    plot1a = scatter!(data_plot_condition[:,:DaysSinceSeeding], data_plot_condition[:,:NecroticRadius],xlim=(0,11), ylim=(0,350), xlab=L"t \mathrm{ \ [days]}",ylab=L"R(t) \ [\mu \mathrm{m}]",xticks = [0,2,4,6,8,10],yticks=[0,100,200,300],titlefont=fnt, guidefont=fnt, tickfont=fnt,markersize = 5, lw=3,legend = false,color=:Black,mode="markers",markerstrokecolor=:Black)
    display(plot1a)
    savefig(plot1a,filepath_save[1] * "Fig1a" * "_" * replace(plot_Conditions[j],":" => "" ) * ".pdf")
end


#######################################################################################
### Plot data - spheroid structure
# Necrotic Radius/Radius (y) v Radius (x)

plot_Conditions =["M100:F0", "M75:F25", "M50:F50","M25:F75"]
for j=1:4
    # filter data for condition j
    data_plot_condition = copy(data_all_v1);
    function filter_plot_condition(Condition)::Bool
        interesting = Condition == plot_Conditions[j]
    end
    filter!([:Condition] => filter_plot_condition, data_plot_condition)
    # plot data and save for each condition
    if j== 1
        global plot1b = scatter(data_plot_condition[:,:Radius],data_plot_condition[:,:NecroticRadius]./data_plot_condition[:,:Radius],xlim=(0,350), ylim=(0,1), xlab=L"R(t) \ [\mu \mathrm{m}]",ylab=L"R_{\mathrm{n}}(t)/R(t) \ [-]",xticks=[0,100,200,300],yticks=[0,0.2,0.4,0.6,0.8,1.0],titlefont=fnt, guidefont=fnt, tickfont=fnt,markersize = 5, lw=3,label=plot_Conditions[j],color=cur_colors[j],mode="markers",markerstrokecolor=cur_colors[j])
    else
        plot1b = scatter!(data_plot_condition[:,:Radius],data_plot_condition[:,:NecroticRadius]./data_plot_condition[:,:Radius],xlim=(0,350), ylim=(0,1), xlab=L"R(t) \ [\mu \mathrm{m}]",ylab=L"R_{\mathrm{n}}(t)/R(t) \ [-]",xticks=[0,100,200,300],yticks=[0,0.2,0.4,0.6,0.8,1.0],titlefont=fnt, guidefont=fnt, tickfont=fnt,markersize = 5, lw=3,label=plot_Conditions[j],color=cur_colors[j],mode="markers",markerstrokecolor=cur_colors[j])
    end
end
display(plot1b)
savefig(plot1b,filepath_save[1] * "Fig1b" * ".pdf")

#######################################################################################
### Loop through each condition for MLE, profiles, and bounds of confidence interval


# There are 4 parameters to estimate R0, s, λ, ℝ

    maxtime_optimisation_MLE = 60; # maximimum time to search for MLE in optimisation
    maxtime_optimisation_profile = 60; # maximum time to search for optima in each step of profiling

    plot_Conditions =["M100:F0", "M75:F25", "M50:F50","M25:F75","M0:F100"]
    for j=1:4
        # filter data for condition j
        data_plot_condition = copy(data_all_v1);
        function filter_plot_condition(Condition)::Bool
            interesting = Condition == plot_Conditions[j]
        end
        filter!([:Condition] => filter_plot_condition, data_plot_condition)
        global t1 = data_plot_condition[:,:DaysSinceSeeding] -2*ones(length(data_plot_condition[:,:DaysSinceSeeding]),1) # measurement times
        global data = hcat(data_plot_condition[:,:Radius], data_plot_condition[:,:NecroticRadius]) # radial measurements

        # initial guesses for each condition
        a=zeros(5)
        if j==1
            RR0=201.0;
            SS=0.24;
            λλ=2.0;
            ℝℝ=225.0;
            Sd=8.0;
        elseif j ==2
            RR0=177.0;
            SS=0.24;
            λλ=2.4;
            ℝℝ=240.0;
            Sd=16.0;
        elseif j == 3
            RR0=177.0;
            SS=0.20;
            λλ=3.0;
            ℝℝ=235.0;
            Sd=17.0;
        elseif j == 4
            RR0=159.0;
            SS=0.20;
            λλ=1.9;
            ℝℝ=210.0;
            Sd=10.0;
        end

        TH=-1.921; #95% confidence interval threshold
        
        function reducedGreenspaneval(t1, a)
            # reduced Greenspan model
            R0=a[1];
            S=a[2]
            λ=a[3];
            ℝ=a[4];
            
            function reducedGreenspan(du,u,p,t)
                R,Rn = u;
                s, λ, ℝ  = p;
                du[1] = (s/3)*(R - Rn^3/R^2) - s*λ*Rn^3/R^2;
                if R > ℝ
                    fRn(x) = R^2 - x^2 - 2*x^2*(1 - x/R) - ℝ^2; # where x=Rn
                    Rnestimate = find_zero(fRn, (0,R)) ;
                    du[2] = Rn - Rnestimate;
                else
                    du[2] = Rn - 0;
                end
                nothing
            end
        M = [1. 0. 
            0.  0.] # mass-matrix as we solve as differential-algebraic system

        f = ODEFunction(reducedGreenspan,mass_matrix=M)
        prob_mm = ODEProblem(f,[R0,0],(0.0,maximum(t1)),(S,λ,ℝ));
        sol = solve(prob_mm,Rodas5(autodiff=false),reltol=1e-8,abstol=1e-8,saveat=t1);
        return(sol.u)

        end

        function error(data,a)
            # error model
            
            # return y values for unique t1
            sort_unique_t1 = sort(unique(t1)); 
            y_tmp=reducedGreenspaneval(sort_unique_t1,a);

            # map unique values of t to t1 vector
            y=zeros(length(t1),2);
            for i=1:length(t1)
                # find t1(i) in sort(unique(t1))
                index_to_lookup = findall(==(t1[i]),sort_unique_t1);
                y[i,:] = y_tmp[index_to_lookup,:][1]
            end
            y_long_format = [y[:,1];y[:,2]];
            data_long_format = [data[:,1];data[:,2]];

            e=0;
            dist=Normal(0,a[5]);
            e=loglikelihood(dist,data_long_format-y_long_format) 
            ee=sum(e)
            return ee
        end

        function fun(a)
            return error(data,a)
        end
        
        # search for MLE for maxtime_optimisation_MLE seconds
        function optimise_MLE(fun,θ₀,lb,ub) 
            tomax = (θ,∂θ) -> fun(θ)
            opt = Opt(:LN_NELDERMEAD,length(θ₀))
            opt.max_objective = tomax
            opt.lower_bounds = lb       # Lower bound
            opt.upper_bounds = ub       # Upper bound
            opt.ftol_abs = 1e-22;
            opt.xtol_abs = 0.0;
            opt.maxtime = maxtime_optimisation_MLE; # maximum time in seconds
            res = optimize(opt,θ₀)
            return res[[2,1]]
        end


        # search for optima for profile  for maxtime_optimisation_profile seconds
        function optimise(fun,θ₀,lb,ub) 
            tomax = (θ,∂θ) -> fun(θ)
            opt = Opt(:LN_NELDERMEAD,length(θ₀))
            opt.max_objective = tomax
            opt.lower_bounds = lb       # Lower bound
            opt.upper_bounds = ub       # Upper bound
            opt.ftol_abs = 1e-22;
            opt.xtol_abs = 0.0;
            opt.maxtime = maxtime_optimisation_profile; # maximum time in seconds
            res = optimize(opt,θ₀)
            return res[[2,1]]
        end
        

        #######################################################################################
        # MLE

        θG = [RR0, SS, λλ, ℝℝ, Sd] # initial guess 
        lb = [150.0, 0.15, 0.001, 190.0, 0.0001] # lower bound
        ub = [210.0, 0.3, 5.0, 260.0, 25.0] #$ upper bound

        # MLE optimisation
        (xopt,fopt)  = optimise_MLE(fun,θG,lb,ub)

        # storing the MLE
        global fmle=fopt
        global R0mle=xopt[1]
        global Smle=xopt[2]
        global λmle=xopt[3]
        global ℝmle=xopt[4]
        global Sdmle=xopt[5]

        # plot the model simulated at the MLE
        t1_smooth = LinRange(minimum(t1),11,101)
        ymle_smooth = reducedGreenspaneval(t1_smooth,xopt);
        ymle_Radius = zeros(length(ymle_smooth),1);
        ymle_NecroticRadius = zeros(length(ymle_smooth),1);
        for i1=1:length(ymle_smooth)
            ymle_Radius[i1] = ymle_smooth[i1,1][1];
            ymle_NecroticRadius[i1] = ymle_smooth[i1,1][2];
        end
        p1 = plot(t1_smooth + 2*ones(length(t1_smooth),1),ymle_Radius,legend=false,xlim=(0,11), ylim=(0,350), xlab=L"t \mathrm{ \ [days]}",ylab=L"R(t) \ [\mu \mathrm{m}]",xticks = [0,2,4,6,8,10],yticks=[0,100,200,300],titlefont=fnt, guidefont=fnt, tickfont=fnt,lw=3,linecolor=:green2)
        p1 = plot!(t1_smooth + 2*ones(length(t1_smooth),1),ymle_NecroticRadius,legend=false,xlim=(0,11), ylim=(0,350), xlab=L"t \mathrm{ \ [days]}",ylab=L"R(t) \ [\mu \mathrm{m}]",xticks = [0,2,4,6,8,10],yticks=[0,100,200,300],titlefont=fnt, guidefont=fnt, tickfont=fnt,lw=3,linecolor=:Black)
        p1 = scatter!(data_plot_condition[:,:DaysSinceSeeding], data_plot_condition[:,:Radius],xlim=(0,11), ylim=(0,350), xlab=L"t \mathrm{ \ [days]}",ylab=L"R(t) \ [\mu \mathrm{m}]",xticks = [0,2,4,6,8,10],yticks=[0,100,200,300],titlefont=fnt, guidefont=fnt, tickfont=fnt,markersize = 5, lw=3,legend = false,color=:green2,mode="markers",markerstrokecolor=:green2)
        p1 = scatter!(data_plot_condition[:,:DaysSinceSeeding], data_plot_condition[:,:NecroticRadius],xlim=(0,11), ylim=(0,350), xlab=L"t \mathrm{ \ [days]}",ylab=L"R(t) \ [\mu \mathrm{m}]",xticks = [0,2,4,6,8,10],yticks=[0,100,200,300],titlefont=fnt, guidefont=fnt, tickfont=fnt,markersize = 5, lw=3,legend = false,color=:Black,mode="markers",markerstrokecolor=:Black) 
        display(p1)
        savefig(p1,filepath_save[1] * "Figp1_" * replace(plot_Conditions[j],":" => "" )  * ".pdf")

        #######################################################################################
        # Profiling (using the MLE as the first guess at each point)
       
        nptss=20;

        #Profile RR0

        R0min=lb[1]
        R0max=ub[1]
        R0range_lower=reverse(LinRange(R0min,R0mle,nptss))
        R0range_upper=LinRange(R0mle + (R0max-R0mle)/nptss,R0max,nptss)
        
        nrange_lower=zeros(4,nptss)
        llR0_lower=zeros(nptss)
        nllR0_lower=zeros(nptss)
        predict_R0_lower=zeros(length(t1_smooth),nptss)
        
        nrange_upper=zeros(4,nptss)
        llR0_upper=zeros(nptss)
        nllR0_upper=zeros(nptss)
        predict_R0_upper=zeros(length(t1_smooth),nptss)
        
        # start at mle and increase parameter (upper)
        for i in 1:nptss
            function fun1(aa)
                return error(data,[R0range_upper[i],aa[1],aa[2],aa[3],aa[4]])
            end
            lb1=[lb[2],lb[3],lb[4],lb[5]];
            ub1=[ub[2],ub[3],ub[4],ub[5]];
            local θG1=[Smle, λmle,ℝmle,Sdmle]
            local (xo,fo)=optimise(fun1,θG1,lb1,ub1)
            nrange_upper[:,i]=xo[:]
            llR0_upper[i]=fo[1]
            if fo > fmle
                println("new MLE R0 upper")
                global fmle=fo
                global R0mle=R0range_upper[i]
                global Smle=nrange_upper[1,i]
                global λmle=nrange_upper[2,i]
                global ℝmle=nrange_upper[3,i]
                global Sdmle=nrange_upper[4,i]
            end
        end
        
        # start at mle and decrease parameter (lower)
        for i in 1:nptss
            function fun1a(aa)
                return error(data,[R0range_lower[i],aa[1],aa[2],aa[3],aa[4]])
            end
            lb1=[lb[2],lb[3],lb[4],lb[5]];
            ub1=[ub[2],ub[3],ub[4],ub[5]];
            local θG1=[Smle, λmle,ℝmle,Sdmle]
            local (xo,fo)=optimise(fun1a,θG1,lb1,ub1)
            nrange_lower[:,i]=xo[:]
            llR0_lower[i]=fo[1]
            if fo > fmle
                println("new MLE R0 lower")
                global fmle = fo
                global R0mle=R0range_lower[i]
                global Smle=nrange_lower[1,i]
                global λmle=nrange_lower[2,i]
                global ℝmle=nrange_lower[3,i]
                global Sdmle=nrange_lower[4,i]
            end
        end
        
        # combine the lower and upper
        R0range = [reverse(R0range_lower); R0range_upper]
        nrange = [reverse(nrange_lower); nrange_upper ]
        llR0 = [reverse(llR0_lower); llR0_upper] 
        nllR0=llR0.-maximum(llR0);
        

        #Profile S
        Smin=lb[2]
        Smax=ub[2]
        Srange_lower=reverse(LinRange(Smin,Smle,nptss))
        Srange_upper=LinRange(Smle + (Smax-Smle)/nptss,Smax,nptss)
        
        nrange_lower=zeros(4,nptss)
        llS_lower=zeros(nptss)
        nllS_lower=zeros(nptss)
        predict_S_lower=zeros(length(t1_smooth),nptss)
        
        nrange_upper=zeros(4,nptss)
        llS_upper=zeros(nptss)
        nllS_upper=zeros(nptss)
        predict_S_upper=zeros(length(t1_smooth),nptss)
        
        # start at mle and increase parameter (upper)
        for i in 1:nptss
            function fun2(aa)
                return error(data,[aa[1],Srange_upper[i],aa[2],aa[3],aa[4]])
            end
            lb1=[lb[1],lb[3],lb[4],lb[5]];
            ub1=[ub[1],ub[3],ub[4],ub[5]];
            local θG1=[R0mle, λmle,ℝmle,Sdmle]
            local (xo,fo)=optimise(fun2,θG1,lb1,ub1)
            nrange_upper[:,i]=xo[:]
            llS_upper[i]=fo[1]
          
            if fo > fmle
                println("new MLE S upper")
                global fmle = fo
                global R0mle=nrange_upper[1,i]
                global Smle=Srange_upper[i]
                global λmle=nrange_upper[2,i]
                global ℝmle=nrange_upper[3,i]
                global Sdmle=nrange_upper[4,i]
            end
        end
        
        # start at mle and decrease parameter (lower)
        for i in 1:nptss
            function fun2a(aa)
                return error(data,[aa[1],Srange_lower[i],aa[2],aa[3],aa[4]])
            end
            lb1=[lb[1],lb[3],lb[4],lb[5]];
            ub1=[ub[1],ub[3],ub[4],ub[5]];
            local θG1=[R0mle, λmle,ℝmle,Sdmle]
            local (xo,fo)=optimise(fun2a,θG1,lb1,ub1)
            nrange_lower[:,i]=xo[:]
            llS_lower[i]=fo[1]
            if fo > fmle
                println("new MLE S lower")
                global fmle = fo
                global R0mle=nrange_lower[1,i]
                global Smle=Srange_lower[i]
                global λmle=nrange_lower[2,i]
                global ℝmle=nrange_lower[3,i]
                global Sdmle=nrange_lower[4,i]
            end
        end
        
        # combine the lower and upper
        Srange = [reverse(Srange_lower);Srange_upper]
        nrange = [reverse(nrange_lower); nrange_upper ]
        llS = [reverse(llS_lower); llS_upper]     
        nllS=llS.-maximum(llS);
        

        #Profile λ
        λmin=lb[3]
        λmax=ub[3]
        λrange_lower=reverse(LinRange(λmin,λmle,nptss))
        λrange_upper=LinRange(λmle + (λmax-λmle)/nptss,λmax,nptss)
        
        nrange_lower=zeros(4,nptss)
        llλ_lower=zeros(nptss)
        nllλ_lower=zeros(nptss)
        predict_λ_lower=zeros(length(t1_smooth),nptss)
        
        nrange_upper=zeros(4,nptss)
        llλ_upper=zeros(nptss)
        nllλ_upper=zeros(nptss)
        predict_λ_upper=zeros(length(t1_smooth),nptss)
        
        # start at mle and increase parameter (upper)
        for i in 1:nptss
            function fun3(aa)
                return error(data,[aa[1],aa[2],λrange_upper[i],aa[3],aa[4]])
            end
            lb1=[lb[1],lb[2],lb[4],lb[5]];
            ub1=[ub[1],ub[2],ub[4],ub[5]];
            local θG1=[R0mle,Smle,ℝmle,Sdmle]    
            local (xo,fo)=optimise(fun3,θG1,lb1,ub1)
            nrange_upper[:,i]=xo[:]
            llλ_upper[i]=fo[1]
            if fo > fmle
                println("new MLE λ upper")
                global fmle = fo
                global R0mle=nrange_upper[1,i]
                global Smle=nrange_upper[2,i]
                global λmle=λrange_upper[i]
                global ℝmle=nrange_upper[3,i]
                global Sdmle=nrange_upper[4,i]
            end
        end
        
        # start at mle and decrease parameter (lower)
        for i in 1:nptss
            function fun3a(aa)
                return error(data,[aa[1],aa[2],λrange_lower[i],aa[3],aa[4]])
            end
            lb1=[lb[1],lb[2],lb[4],lb[5]];
            ub1=[ub[1],ub[2],ub[4],ub[5]];
            local θG1=[R0mle,Smle,ℝmle,Sdmle]      
            local (xo,fo)=optimise(fun3a,θG1,lb1,ub1)
            nrange_lower[:,i]=xo[:]
            llλ_lower[i]=fo[1]
            if fo > fmle
                println("new MLE λ lower")
                global fmle = fo
                global R0mle=nrange_lower[1,i]
                global Smle=nrange_lower[2,i]
                global λmle=λrange_lower[i]
                global ℝmle=nrange_lower[3,i]
                global Sdmle=nrange_lower[4,i]
            end
        end
        
        # combine the lower and upper
        λrange = [reverse(λrange_lower);λrange_upper]
        nrange = [reverse(nrange_lower); nrange_upper ]
        llλ = [reverse(llλ_lower); llλ_upper] 
        
        nllλ=llλ.-maximum(llλ)

        
        #Profile ℝ
        ℝmin=lb[4]
        ℝmax=ub[4]
        ℝrange_lower=reverse(LinRange(ℝmin,ℝmle,nptss))
        ℝrange_upper=LinRange(ℝmle + (ℝmax-ℝmle)/nptss,ℝmax,nptss)
        
        nrange_lower=zeros(4,nptss)
        llℝ_lower=zeros(nptss)
        nllℝ_lower=zeros(nptss)
        predict_ℝ_lower=zeros(length(t1_smooth),nptss)
        
        nrange_upper=zeros(4,nptss)
        llℝ_upper=zeros(nptss)
        nllℝ_upper=zeros(nptss)
        predict_ℝ_upper=zeros(length(t1_smooth),nptss)
        
        # start at mle and increase parameter (upper)
        for i in 1:nptss
            function fun4(aa)
                return error(data,[aa[1],aa[2],aa[3],ℝrange_upper[i],aa[4]])
            end
            lb1=[lb[1],lb[2],lb[3],lb[5]];
            ub1=[ub[1],ub[2],ub[3],ub[5]];
            local θG1=[R0mle,Smle,λmle,Sdmle]
            local (xo,fo)=optimise(fun4,θG1,lb1,ub1)
            nrange_upper[:,i]=xo[:]
            llℝ_upper[i]=fo[1]
            if fo > fmle
                println("new MLE ℝ upper")
                global fmle = fo
                global R0mle=nrange_upper[1,i]
                global Smle=nrange_upper[2,i]
                global λmle=nrange_upper[3,i]
                global ℝmle=ℝrange_upper[i]
                global Sdmle=nrange_upper[4,i]
            end
        end
        
        # start at mle and decrease parameter (lower)
        for i in 1:nptss
            function fun4a(aa)
                return error(data,[aa[1],aa[2],aa[3],ℝrange_lower[i],aa[4]])
            end
            lb1=[lb[1],lb[2],lb[3],lb[5]];
            ub1=[ub[1],ub[2],ub[3],ub[5]];
            local θG1=[R0mle,Smle,λmle,Sdmle]
            local (xo,fo)=optimise(fun4a,θG1,lb1,ub1)
            nrange_lower[:,i]=xo[:]
            llℝ_lower[i]=fo[1]
            if fo > fmle
                println("new MLE ℝ lower")
                global fmle = fo
                global R0mle=nrange_lower[1,i]
                global Smle=nrange_lower[2,i]
                global λmle=nrange_lower[3,i]
                global ℝmle=ℝrange_lower[i]
                global Sdmle=nrange_lower[4,i]
            end
        end
        
        # combine the lower and upper
        ℝrange = [reverse(ℝrange_lower);ℝrange_upper]
        llℝ = [reverse(llℝ_lower); llℝ_upper]        
        nllℝ=llℝ.-maximum(llℝ);
        
        
        #Profile Sd
        Sdmin=lb[5]+0.001
        Sdmax=ub[5]
        Sdrange_lower=reverse(LinRange(Sdmin,Sdmle,nptss))
        Sdrange_upper=LinRange(Sdmle + (Sdmax-Sdmle)/nptss,Sdmax,nptss)
        
        nrange_lower=zeros(4,nptss)
        llSd_lower=zeros(nptss)
        nllSd_lower=zeros(nptss)
        predict_Sd_lower=zeros(length(t1_smooth),nptss)
        
        nrange_upper=zeros(4,nptss)
        llSd_upper=zeros(nptss)
        nllSd_upper=zeros(nptss)
        predict_Sd_upper=zeros(length(t1_smooth),nptss)
        
        # start at mle and increase parameter (upper)
        for i in 1:nptss
            function fun5(aa)
                return error(data,[aa[1],aa[2],aa[3],aa[4],Sdrange_upper[i],])
            end
            lb1=[lb[1],lb[2],lb[3],lb[4]];
            ub1=[ub[1],ub[2],ub[3],ub[4]];
            local θG1=[R0mle,Smle,λmle,ℝmle]
            local (xo,fo)=optimise(fun5,θG1,lb1,ub1)
            nrange_upper[:,i]=xo[:]
            llSd_upper[i]=fo[1]
            if fo > fmle
                println("new MLE Sd upper")
                global fmle = fo
                global R0mle=nrange_upper[1,i]
                global Smle=nrange_upper[2,i]
                global λmle=nrange_upper[3,i]
                global ℝmle=nrange_upper[4,i]
                global Sdmle=Sdrange_upper[i]
            end
        end
        
        # start at mle and decrease parameter (lower)
        for i in 1:nptss
            function fun5a(aa)
                return error(data,[aa[1],aa[2],aa[3],aa[4],Sdrange_lower[i],])
            end
            lb1=[lb[1],lb[2],lb[3],lb[4]];
            ub1=[ub[1],ub[2],ub[3],ub[4]];
            local θG1=[R0mle,Smle,λmle,ℝmle]
            local (xo,fo)=optimise(fun5a,θG1,lb1,ub1)
            nrange_lower[:,i]=xo[:]
            llSd_lower[i]=fo[1]
            if fo > fmle
                println("new MLE Sd lower")
                global fmle = fo
                global R0mle=nrange_upper[1,i]
                global Smle=nrange_upper[2,i]
                global λmle=nrange_upper[3,i]
                global ℝmle=nrange_upper[4,i]
                global Sdmle=Sdrange_upper[i]
            end
        end
        
        # combine the lower and upper
        Sdrange = [reverse(Sdrange_lower);Sdrange_upper]
        nrange = [reverse(nrange_lower); nrange_upper ]
        llSd = [reverse(llSd_lower); llSd_upper] 
        nllSd=llSd.-maximum(llSd);

        ##############################################################################################################
        combined_plot_linewidth = 2;

        # interpolate for smoother profile likelihoods
        interp_nptss= 1001;
        
        # R0
        interp_points_R0range =  LinRange(R0min,R0max,interp_nptss)
        interp_R0 = LinearInterpolation(R0range,nllR0)
        interp_nllR0 = interp_R0(interp_points_R0range)
        
        s1=plot(interp_points_R0range,interp_nllR0,xlim=(R0min,R0max),ylim=(-4,0.1),yticks=[-3,-2,-1,0],xlab=L"R(2) \ [\mu \mathrm{m}]",ylab=L"\hat{\ell}_{p}",legend=false,lw=5,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=:deepskyblue3)
        s1=hline!([-1.92],lw=2,linecolor=:black,linestyle=:dot)
        s1=vline!([R0mle],lw=3,linecolor=:red)

        if j==1
        global s1_combined=hline([-1.92],lw=2,linecolor=:black,linestyle=:dot)
        end
        s1_combined=plot!(s1_combined,interp_points_R0range,interp_nllR0,xlim=(R0min,R0max),ylim=(-4,0.1),yticks=[-3,-2,-1,0],xlab=L"R(2) \  [\mu \mathrm{m}]",ylab=L"\hat{\ell}_{p}",legend=false,lw=combined_plot_linewidth,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=cur_colors[j])
        s1_combined=vline!(s1_combined,[R0mle],lw=combined_plot_linewidth,linecolor=cur_colors[j],linestyle=:dash)

        
        # S
        interp_points_Srange =  LinRange(Smin,Smax,interp_nptss)
        interp_S = LinearInterpolation(Srange,nllS)
        interp_nllS = interp_S(interp_points_Srange)
        
        s2=plot(interp_points_Srange,interp_nllS,xlim=(Smin,Smax),ylim=(-4,0.1),yticks=[-3,-2,-1,0],xlab=L"s \ [\mathrm{day}^{-1}]",ylab=L"\hat{\ell}_{p}",legend=false,lw=5,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=:deepskyblue3)
        s2=hline!([-1.92],lw=2,linecolor=:black,linestyle=:dot)
        s2=vline!([Smle],lw=3,linecolor=:red)

        if j==1
        global  s2_combined=hline([-1.92],lw=2,linecolor=:black,linestyle=:dot)
        end
        s2_combined=plot!(s2_combined,interp_points_Srange,interp_nllS,xlim=(Smin,Smax),ylim=(-4,0.1),yticks=[-3,-2,-1,0],xlab=L"s \ [\mathrm{day}^{-1}]",ylab=L"\hat{\ell}_{p}",legend=false,lw=combined_plot_linewidth,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=cur_colors[j])
        s2_combined=vline!(s2_combined,[Smle],lw=combined_plot_linewidth,linecolor=cur_colors[j],linestyle=:dash)

        
        # λ
        
        interp_points_λrange =  LinRange(λmin,λmax,interp_nptss)
        interp_λ = LinearInterpolation(λrange,nllλ)
        interp_nllλ = interp_λ(interp_points_λrange)
        
        s3=plot(interp_points_λrange,interp_nllλ,xlim=(λmin,λmax),ylim=(-4,0.1),yticks=[-3,-2,-1,0],xlab=L"\gamma \ [-]",ylab=L"\hat{\ell}_{p}",legend=false,lw=5,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=:deepskyblue3)
        s3=hline!([-1.92],lw=2,linecolor=:black,linestyle=:dot)
        s3=vline!([λmle],lw=3,linecolor=:red)

        if j==1
        global  s3_combined=hline([-1.92],lw=2,linecolor=:black,linestyle=:dot)
        end
        s3_combined=plot!(s3_combined,interp_points_λrange,interp_nllλ,xlim=(λmin,λmax),ylim=(-4,0.1),yticks=[-3,-2,-1,0],xlab=L"\gamma \ [-]",ylab=L"\hat{\ell}_{p}",legend=false,lw=combined_plot_linewidth,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=cur_colors[j])
            s3_combined=vline!(s3_combined,[λmle],lw=combined_plot_linewidth,linecolor=cur_colors[j],linestyle=:dash)
        

        # ℝ
        
        interp_points_ℝrange =  LinRange(ℝmin,ℝmax,interp_nptss)
        interp_ℝ = LinearInterpolation(ℝrange,nllℝ)
        interp_nllℝ = interp_ℝ(interp_points_ℝrange)
        
        s4=plot(interp_points_ℝrange,interp_nllℝ,xlim=(ℝmin,ℝmax),ylim=(-4,0.1),yticks=[-3,-2,-1,0],xlab=L"\mathcal{R} \ [\mu \mathrm{m}]",ylab=L"\hat{\ell}_{p}",legend=false,lw=5,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=:deepskyblue3)
        s4=hline!([-1.92],lw=2,linecolor=:black,linestyle=:dot)
        s4=vline!([ℝmle],lw=3,linecolor=:red)

        if j==1
        global  s4_combined=hline([-1.92],lw=2,linecolor=:black,linestyle=:dot)
        end
        s4_combined=plot!(s4_combined,interp_points_ℝrange,interp_nllℝ,xlim=(ℝmin,ℝmax),ylim=(-4,0.1),yticks=[-3,-2,-1,0],xlab=L"\mathcal{R} \ [\mu \mathrm{m}]",ylab=L"\hat{\ell}_{p}",legend=false,lw=combined_plot_linewidth,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=cur_colors[j])
        s4_combined=vline!(s4_combined,[ℝmle],lw=combined_plot_linewidth,linecolor=cur_colors[j],linestyle=:dash)

        
        # Sd
        
        interp_points_Sdrange =  LinRange(Sdmin,Sdmax,interp_nptss)
        interp_Sd = LinearInterpolation(Sdrange,nllSd)
        interp_nllSd = interp_Sd(interp_points_Sdrange)
        
        s5=plot(interp_points_Sdrange,interp_nllSd,xlim=(Sdmin,Sdmax),ylim=(-4,0.1),yticks=[-3,-2,-1,0],xlab=L"\sigma [\mu \mathrm{m}]",ylab=L"\hat{\ell}_{p}",legend=false,lw=5,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=:deepskyblue3)
        s5=hline!([-1.92],lw=2,linecolor=:black,linestyle=:dot)
        s5=vline!([Sdmle],lw=3,linecolor=:red)
        
        if j==1
            global  s5_combined=hline([-1.92],lw=2,linecolor=:black,linestyle=:dot)
        end
        s5_combined=plot!(s5_combined,interp_points_Sdrange,interp_nllSd,xlim=(Sdmin,Sdmax),ylim=(-4,0.1),yticks=[-3,-2,-1,0],xlab=L"\sigma \  [\mu \mathrm{m}]",ylab=L"\hat{\ell}_{p}",legend=false,lw=combined_plot_linewidth,titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt,linecolor=cur_colors[j])
        s5_combined=vline!(s5_combined,[Sdmle],lw=combined_plot_linewidth,linecolor=cur_colors[j],linestyle=:dash)

        # Recompute best fit
        t1_smooth = LinRange(minimum(t1),11,101)
        ymle_smooth = reducedGreenspaneval(t1_smooth,[R0mle;Smle;λmle;ℝmle;Sdmle]);
        ymle_Radius = zeros(length(ymle_smooth),1);
        ymle_NecroticRadius = zeros(length(ymle_smooth),1);
        for i1=1:length(ymle_smooth)
            ymle_Radius[i1] = ymle_smooth[i1,1][1];
            ymle_NecroticRadius[i1] = ymle_smooth[i1,1][2];
        end
        p1_updated = plot(t1_smooth + 2*ones(length(t1_smooth),1),ymle_Radius,legend=false,xlim=(0,11), ylim=(0,350), xlab=L"t \mathrm{ \ [days]}",ylab=L"R(t) \ [\mu \mathrm{m}]",xticks = [0,2,4,6,8,10],yticks=[0,100,200,300],titlefont=fnt, guidefont=fnt, tickfont=fnt,lw=3,linecolor=:green2)
        p1_updated = plot!(t1_smooth + 2*ones(length(t1_smooth),1),ymle_NecroticRadius,legend=false,xlim=(0,11), ylim=(0,350), xlab=L"t \mathrm{ \ [days]}",ylab=L"R(t) \ [\mu \mathrm{m}]",xticks = [0,2,4,6,8,10],yticks=[0,100,200,300],titlefont=fnt, guidefont=fnt, tickfont=fnt,lw=3,linecolor=:Black)
        p1_updated = scatter!(data_plot_condition[:,:DaysSinceSeeding], data_plot_condition[:,:Radius],xlim=(0,11), ylim=(0,350), xlab=L"t \mathrm{ \ [days]}",ylab=L"R(t) \ [\mu \mathrm{m}]",xticks = [0,2,4,6,8,10],yticks=[0,100,200,300],titlefont=fnt, guidefont=fnt, tickfont=fnt,markersize = 5, lw=3,legend = false,color=:green2,mode="markers",markerstrokecolor=:green2)
        p1_updated = scatter!(data_plot_condition[:,:DaysSinceSeeding], data_plot_condition[:,:NecroticRadius],xlim=(0,11), ylim=(0,350), xlab=L"t \mathrm{ \ [days]}",ylab=L"R(t) \ [\mu \mathrm{m}]",xticks = [0,2,4,6,8,10],yticks=[0,100,200,300],titlefont=fnt, guidefont=fnt, tickfont=fnt,markersize = 5, lw=3,legend = false,color=:Black,mode="markers",markerstrokecolor=:Black) 

        # save figures
        display(p1_updated)
        savefig(p1_updated,filepath_save[1] * "Figp1_updated" * replace(plot_Conditions[j],":" => "" )   * ".pdf")
        
        display(s1)
        savefig(s1,filepath_save[1] * "Figs1" * replace(plot_Conditions[j],":" => "" )   * ".pdf")
        
        display(s2)
        savefig(s2,filepath_save[1] * "Figs2" * replace(plot_Conditions[j],":" => "" )   * ".pdf")
        
        display(s3)
        savefig(s3,filepath_save[1] * "Figs3" * replace(plot_Conditions[j],":" => "" )   * ".pdf")
        
        display(s4)
        savefig(s4,filepath_save[1] * "Figs4" * replace(plot_Conditions[j],":" => "" )   * ".pdf")
        
        display(s5)
        savefig(s5,filepath_save[1] * "Figs5" * replace(plot_Conditions[j],":" => "" )   * ".pdf")
        
        display(s1_combined)
        savefig(s1_combined,filepath_save[1] * "Figs1_combined" * replace(plot_Conditions[j],":" => "" )   * ".pdf")
        
        display(s2_combined)
        savefig(s2_combined,filepath_save[1] * "Figs2_combined" * replace(plot_Conditions[j],":" => "" )   * ".pdf")
        
        display(s3_combined)
        savefig(s3_combined,filepath_save[1] * "Figs3_combined" * replace(plot_Conditions[j],":" => "" )   * ".pdf")
        
        display(s4_combined)
        savefig(s4_combined,filepath_save[1] * "Figs4_combined" * replace(plot_Conditions[j],":" => "" )   * ".pdf")
        
        display(s5_combined)
        savefig(s5_combined,filepath_save[1] * "Figs5_combined" * replace(plot_Conditions[j],":" => "" )   * ".pdf")



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
        
        # R0
        (lb_CI_R0,ub_CI_R0) = fun_interpCI(R0mle,interp_points_R0range,interp_nllR0,TH)
        println(round(lb_CI_R0; digits = 4))
        println(round(ub_CI_R0; digits = 4))
        
        # S
        (lb_CI_S,ub_CI_S) = fun_interpCI(Smle,interp_points_Srange,interp_nllS,TH)
        println(round(lb_CI_S; digits = 3))
        println(round(ub_CI_S; digits = 3))
        
        # λ
        (lb_CI_λ,ub_CI_λ) = fun_interpCI(λmle,interp_points_λrange,interp_nllλ,TH)
        println(round(lb_CI_λ; digits = 3))
        println(round(ub_CI_λ; digits = 3))
        
        # ℝ
        (lb_CI_ℝ,ub_CI_ℝ) = fun_interpCI(ℝmle,interp_points_ℝrange,interp_nllℝ,TH)
        println(round(lb_CI_ℝ; digits = 3))
        println(round(ub_CI_ℝ; digits = 3))
        
        # Sd
        (lb_CI_Sd,ub_CI_Sd) = fun_interpCI(Sdmle,interp_points_Sdrange,interp_nllSd,TH)
        println(round(lb_CI_Sd; digits = 3))
        println(round(ub_CI_Sd; digits = 3))


        # Export MLE and bounds to csv (one file for all data) -- -MLE ONLY 
        if @isdefined(df_MLEBoundsAll) == 0
            println("not defined")
            global df_MLEBoundsAll = DataFrame(Condition = replace(plot_Conditions[j],":" => "" ),R0mle=R0mle, lb_CI_R0=lb_CI_R0, ub_CI_R0=ub_CI_R0, Smle=Smle, lb_CI_S=lb_CI_S, ub_CI_S=ub_CI_S, gammamle=λmle, lb_CI_gamma=lb_CI_λ, ub_CI_gamma=ub_CI_λ, bbRmle=ℝmle, lb_CI_bbR=lb_CI_ℝ, ub_CI_bbR=ub_CI_ℝ, Sdmle=Sdmle, lb_CI_Sd=lb_CI_Sd, ub_CI_Sd=ub_CI_Sd )
        else 
            println("DEFINED")
            global df_MLEBoundsAll_thisrow = DataFrame(Condition = replace(plot_Conditions[j],":" => "" ), R0mle=R0mle, lb_CI_R0=lb_CI_R0, ub_CI_R0=ub_CI_R0, Smle=Smle, lb_CI_S=lb_CI_S, ub_CI_S=ub_CI_S, gammamle=λmle, lb_CI_gamma=lb_CI_λ, ub_CI_gamma=ub_CI_λ, bbRmle=ℝmle, lb_CI_bbR=lb_CI_ℝ, ub_CI_bbR=ub_CI_ℝ, Sdmle=Sdmle, lb_CI_Sd=lb_CI_Sd, ub_CI_Sd=ub_CI_Sd)
            append!(df_MLEBoundsAll,df_MLEBoundsAll_thisrow)
        end
        
        CSV.write(filepath_save[1] * "MLEBoundsALL.csv", df_MLEBoundsAll)

end