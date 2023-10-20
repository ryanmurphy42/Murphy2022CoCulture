# Murphy et al. (2022) - CoCulture Spheroids - Plots for Figure 10A and 10B (TWO population TWO compartment reduced Greenspan model with heterogeneous s)

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

pyplot() # plot options
fnt = Plots.font("sans-serif", 20) # plot options
global cur_colors = palette(:default) # plot options
isdir(pwd() * "\\Fig010AB\\") || mkdir(pwd() * "\\Fig010AB") # make folder to save figures if doesnt already exist
filepath_save = [pwd() * "\\Fig010AB\\"] # location to save figures

#######################################################################################
## Load data
data_all_v1 = CSV.read(pwd() * "\\Fig010data.csv", DataFrame);

#######################################################################################
### Plot data - scatter plot 

    plot_Conditions =["M100:F0", "M75:F25", "M50:F50","M25:F75","M0:F100"]
    j=3 # focusing on M50:F50 spheroid data
        
    data_plot_condition = copy(data_all_v1);
    function filter_plot_condition(Condition)::Bool
        interesting = Condition == plot_Conditions[j]
    end
    filter!([:Condition] => filter_plot_condition, data_plot_condition)

    global t1_data = data_plot_condition[:,:DaysSinceSeeding] -2*ones(length(data_plot_condition[:,:DaysSinceSeeding]),1) # measurment times
    global data = hcat(data_plot_condition[:,:Radius], data_plot_condition[:,:NecroticRadius]) # radial measurements

    # plot
    plot_data = scatter(data_plot_condition[:,:DaysSinceSeeding], data_plot_condition[:,:Radius],xlim=(0,11), ylim=(0,350), xlab=L"t \mathrm{ \ [days]}",ylab=L"R(t) \ [\mu \mathrm{m}]",xticks = [0,2,4,6,8,10],yticks=[0,100,200,300],titlefont=fnt, guidefont=fnt, tickfont=fnt,markersize = 5, lw=3,legend = false,color=:green2,mode="markers",markerstrokecolor=:green2)
    plot_data = scatter!(plot_data,data_plot_condition[:,:DaysSinceSeeding], data_plot_condition[:,:NecroticRadius],xlim=(0,11), ylim=(0,350), xlab=L"t \mathrm{ \ [days]}",ylab=L"R(t) \ [\mu \mathrm{m}]",xticks = [0,2,4,6,8,10],yticks=[0,100,200,300],titlefont=fnt, guidefont=fnt, tickfont=fnt,markersize = 5, lw=3,legend = false,color=:Black,mode="markers",markerstrokecolor=:Black) 


############################## MODEL SIMULATION ########################################
# simulate the two compartment two population model
# parameter values are chosen to capture, where possible, loss of fibroblasts from the periphery and chosen so that model simulations agree with measurements of $R(t)$ and $R_{n}(t)$ from the M50:F50 condition 
# therefore, to fit parameters we optimise around an initial guess that we choose base on the one population reduced Greenspan model results 

# For context, see Supplementary Section S5.2.2

    # PHASE I
    function model_phasei(du,u,p,t)
        #Vp1c1, Vp2c1, Vp1c2, Vp2c2, Vp1c3, Vp2c3, R1, R2, R, gc1, gc2 = u # variables
        Rd1, Rd2, s1, s2, gamma_1, gamma_2, omega = p # parameters

            fvp1c1 = (4*pi/3)*s1*(u[9]^3 - u[7]^3)*(u[1]/(u[1]+u[2]));
            fvp2c1 = (4*pi/3)*s2*(u[9]^3 - u[7]^3)*(u[2]/(u[1]+u[2]));
            fvp1c2 = 0;
            fvp2c2 = 0;
            fvp1c3 = 0;
            fvp2c3 = 0;
            g1c1 = 0;
            g1c2 = 0;

            du[1] = fvp1c1 - u[10]*g1c1; # dVp1c1/dt
            du[2] = fvp2c1 - u[10]*(1-g1c1); # dVp2c1/dt
            du[3] = fvp1c2 + u[10]*g1c1 - u[11]*g1c2; # dVp1c2/dt
            du[4] = fvp2c2 + u[10]*(1-g1c1) - u[11]*(1-g1c2); # dVp2c2/dt
            du[5] = fvp1c3 - u[11]*g1c2; # dVp1c3/dt
            du[6] = fvp2c3 - u[11]*(1-g1c2); # dVp2c3/dt
            du[9] = (1/((4*pi/3)*3*u[9]^2))*( fvp1c1 + fvp2c1 + fvp1c2 + fvp2c2 + fvp1c3 + fvp2c3); # dR/dt
            du[10] = -u[10] + 0; #vp1c1 + fvp2c1 - (4*pi/3)*du[9]*(3*u[9]^2  - 3*u[7]^2*dR1dRt ); # g1
            du[11] = -u[11] + 0; #fvp1c2 + fvp2c2 + u[10] - (4*pi/3)*du[9]*(3*u[7]^2 *dR1dRt - 3*u[8]^2*dR2dRt ); # g2

            du[7] = u[7]-0; # R1
            du[8] = u[8]-0; # R2
    end

    # PHASE II
    function model_phaseii(du,u,p,t)
        
        #Vp1c1, Vp2c1, Vp1c2, Vp2c2, Vp1c3, Vp2c3, R1, R2, R, gc1, gc2 = u 
        Rd1, Rd2, s1, s2, gamma_1,gamma_2,omega = p

            h = omega*(4*pi/3)*(u[7]^3 - u[8]^3)*(u[3]/(u[3]+u[4]))*(4*pi/3)*s1*(4*pi/3)*(u[9]^3 - u[7]^3)*(u[2]/(u[1]+u[2]));
  
            fvp1c1 = (4*pi/3)*s1*(u[9]^3 - u[7]^3)*(u[1]/(u[1]+u[2])) + h;
            fvp2c1 = (4*pi/3)*s2*(u[9]^3 - u[7]^3)*(u[2]/(u[1]+u[2])) - h;
            fvp1c2 = -(4*pi/3)*s1*3*gamma_1*(u[7]^3 - u[8]^3)*(u[3]/(u[3]+u[4])) - h;
            fvp2c2 = -(4*pi/3)*s1*3*gamma_2*(u[7]^3 - u[8]^3)*(u[4]/(u[3]+u[4])) - h;
            fvp1c3 = 0;
            fvp2c3 = 0;
            g1c1 = (u[1]/(u[1]+u[2]));
            g1c2 = 0;
          
            dR1dRt = (2*u[9] + 2*(u[7]^2/u[9]^2)*(u[9] - u[7]) - 2*u[7]^2 /u[9] )./(2*u[7] + 4*(u[7]/u[9])*(u[9] - u[7]) - 2*u[7]^2 /u[9]);

            du[1] = fvp1c1 - u[10]*g1c1; # dVp1c1/dt
            du[2] = fvp2c1 - u[10]*(1-g1c1); # dVp2c1/dt
            du[3] = fvp1c2 + u[10]*g1c1 - u[11]*g1c2; # dVp1c2/dt
            du[4] = fvp2c2 + u[10]*(1-g1c1) - u[11]*(1-g1c2); # dVp2c2/dt
            du[5] = fvp1c3 - u[11]*g1c2; # dVp1c3/dt
            du[6] = fvp2c3 - u[11]*(1-g1c2); # dVp2c3/dt
            du[9] = (1/((4*pi/3)*3*u[9]^2))*( fvp1c1 + fvp2c1 + fvp1c2 + fvp2c2 + fvp1c3 + fvp2c3); # dR/dt
            du[10] = -u[10] + fvp1c1 + fvp2c1 - (4*pi/3)*du[9]*(3*u[9]^2 - 3*u[7]^2*dR1dRt ); # g1
            du[11] = -u[11] + 0; # fvp1c2 + fvp2c2 + u[10] - (4*pi/3)*du[9]*(3*u[7]^2 *dR1dRt - 3*u[8]^2*dR2dRt] ); # g2

            du[7] = -Rd1^2 + u[9]^2 - u[7]^2 - 2*(u[7]^2 /u[9])*( u[9] - u[7]); # R1
            du[8] = u[8]-0; # R2
        
    end
    
    # Mass matrix for solving differential-algebraic system

    M = [1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. # Vp1c1 = u[1]
        0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0.# Vp2c1 = u[2]
        0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0.# Vp1c2 = u[3]
        0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.# Vp2c2 = u[4]
        0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0.# Vp1c3 = u[5]
        0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0.# Vp2c3 = u[6]
        0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.# R1 = u[7]
        0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.# R2 = u[8]
        0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0.# R = u[9]
        0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.# g1 = u[10]
        0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.];# g2 = u[11]



function twopoptwocompartmentreducedGreenspan(t1, a)
    # function to simulate the two population two compartment reduced Greenspan model


    # variables   
    # Vp1c1 = u[1]
    # Vp2c1 = u[2]
    # Vp1c2 = u[3]
    # Vp2c2 = u[4]
    # Vp1c3 = u[5]
    # Vp2c3 = u[6]
    # R1 = u[7]
    # R2 = u[8]
    # R = u[9]
    # g1  u[10]
    # g2 = u[11]

    # initial conditions and parameters
    initial_proportion_pop1 = 0.5;
    R_init_phasei = a[1];
    s1 = a[2];
    s2 = a[3];
    gamma_1= a[4]; 
    gamma_2= a[4]; 
    global Rd1 = a[5];
    global Rd2 = 10000; # not used here
    omega = 0.0;

            ###### solve phasei

            # initial conditions
            Vp1c1_init_phasei = initial_proportion_pop1*(4*pi/3)*R_init_phasei^3;
            Vp2c1_init_phasei = (4*pi/3)*R_init_phasei^3 - Vp1c1_init_phasei;
            Vp1c2_init_phasei = 0;
            Vp2c2_init_phasei = 0;
            Vp1c3_init_phasei = 0;
            Vp2c3_init_phasei = 0;
            R1_init_phasei = 0;
            R2_init_phasei = 0;
            g1_init_phasei = 0;
            g2_init_phasei = 0;
            u0model_phasei =  [Vp1c1_init_phasei,Vp2c1_init_phasei, Vp1c2_init_phasei, Vp2c2_init_phasei,Vp1c3_init_phasei,Vp2c3_init_phasei,R1_init_phasei,R2_init_phasei,R_init_phasei,g1_init_phasei,g2_init_phasei];

            p_vec = (Rd1,Rd2,s1, s2,gamma_1, gamma_2,omega) # parameters

            # stopping condition for phase i
            function condition_terminate_phasei(u,t,integrator)
                u[9] - Rd1 -1e-2
            end 
            affect!(integrator) = terminate!(integrator)
            cb_terminate_phasei = ContinuousCallback(condition_terminate_phasei, affect!)


            t1_phasei = LinRange(0.0,maximum(t1),101) # time to simulate model
            tspan_phasei = (minimum(t1_phasei),maximum(t1_phasei)) # time to simulate model
            fmodel_phasei = ODEFunction(model_phasei,mass_matrix=M) # ODEfunction to simulate
            probmodel_phasei = ODEProblem(fmodel_phasei,u0model_phasei,tspan_phasei,p_vec) # ODEproblem to simulate
            solmodel_phasei = solve(probmodel_phasei,Rodas5(autodiff=false),reltol=1e-8,abstol=1e-8,saveat=t1,callback=cb_terminate_phasei); # simulating the model


            ######  solve phase ii
            
            #initial conditions
            Vp1c1_init_phaseii = solmodel_phasei.u[end][1];
            Vp2c1_init_phaseii = solmodel_phasei.u[end][2];
            R_init_phaseii = solmodel_phasei.u[end][9];
            # computation to compute initial condition for R1_init_phaseii
            fRn(x) = R_init_phaseii^2 - x^2 - 2*x^2*(1 - x/R_init_phaseii) - Rd1^2; # where x=Rn
            Rnestimate = find_zero(fRn, (0,R_init_phaseii)) ;
            R1_init_phaseii = Rnestimate;
            #other initial conditions 
            R2_init_phaseii = 0;
            g1_init_phaseii = 0;
            g2_init_phaseii = 0;
            Vp1c2_init_phaseii = (4*pi/3)*R1_init_phaseii^3 *(Vp1c1_init_phaseii/(Vp1c1_init_phaseii + Vp2c1_init_phaseii));
            Vp2c2_init_phaseii = (4*pi/3)*R1_init_phaseii^3 *(Vp2c1_init_phaseii/(Vp1c1_init_phaseii + Vp2c1_init_phaseii));
            Vp1c3_init_phaseii = 0;
            Vp2c3_init_phaseii = 0;
            u0model_phaseii =  [Vp1c1_init_phaseii,Vp2c1_init_phaseii, Vp1c2_init_phaseii, Vp2c2_init_phaseii,Vp1c3_init_phaseii,Vp2c3_init_phaseii,R1_init_phaseii,R2_init_phaseii,R_init_phaseii,g1_init_phaseii,g2_init_phaseii];

            # stopping condition for phase ii (only two compartments so Rd2 chosen so it terminates at max time)
            function condition_terminate_phaseii(u,t,integrator)
                u[9] - Rd2
            end 
            affect!(integrator) = terminate!(integrator)
            cb_terminate_phaseii = ContinuousCallback(condition_terminate_phaseii, affect!)

            t1_phaseii = LinRange(solmodel_phasei.t[end],maximum(t1),101) # time to simulate model
            tspan_phaseii = (minimum(t1_phaseii),maximum(t1_phaseii)) # time to simulate model
            fmodel_phaseii = ODEFunction(model_phaseii,mass_matrix=M) # ODEfunction to simulate
            probmodel_phaseii = ODEProblem(fmodel_phaseii,u0model_phaseii,tspan_phaseii,p_vec)  # ODEproblem to simulate
            solmodel_phaseii = solve(probmodel_phaseii,Rodas5(autodiff=false),reltol=1e-8,abstol=1e-8,saveat=t1,callback=cb_terminate_phaseii); # simulating the model

        return([[solmodel_phasei.t;solmodel_phaseii.t] [solmodel_phasei.u;solmodel_phaseii.u ]])

end



function error(data,a)
    # error model  

    # return y values for unique t1
    sort_unique_t1 = sort(unique(t1_data));
    return_tmp=twopoptwocompartmentreducedGreenspan(sort_unique_t1,a);
    t_tmp = return_tmp[:,1]; # times 
    y_tmp = return_tmp[:,2]; # variables
    # map unique values of t to t1 vector
    y=zeros(length(t1_data),2);
    for i=1:length(t1_data)
        # find t1(i) in t_tmp
        index_to_lookup = findall(==(t1_data[i]),t_tmp);
        y[i,:] = [y_tmp[index_to_lookup][1][9] y_tmp[index_to_lookup][1][7]]
    end
    y_long_format = [y[:,1];y[:,2]];
    data_long_format = [data[:,1];data[:,2]];

    e=0;
    dist=Normal(0,a[6]);
    e=loglikelihood(dist,data_long_format-y_long_format) 
    ee=sum(e)
    return ee
end

function fun(a)
    return error(data,a)
end

# search for optima for 60 seconds (may not be the MLE due to higher-dimensional parameter space)
function optimise_MLE(fun,θ₀,lb,ub) 
    tomax = (θ,∂θ) -> fun(θ)
    opt = Opt(:LN_NELDERMEAD,length(θ₀))
    opt.max_objective = tomax
    opt.lower_bounds = lb       # Lower bound
    opt.upper_bounds = ub       # Upper bound
    opt.ftol_abs = 1e-10;
    opt.xtol_abs = 0.0;
    opt.maxtime = 60; # maximum time in seconds
    res = optimize(opt,θ₀)
    return res[[2,1]]
end

# initial guess for parameters chosen to capture, where possible, loss of fibroblasts from the periphery and chosen so that model simulations agree with measurements of $R(t)$ and $R_{n}(t)$ from the M50:F50 condition 
RR0=177;
SS1=0.6;
SS2=0.0;
λλ=3.07;
ℝℝ=235.0;
Sd=17.0;

θG = [RR0, SS1, SS2, λλ, ℝℝ, Sd] # first guess
lb = [150.0, 0.1, 0.0, 0.001, 180.0, 0.0001] # lower bound
ub = [210.0, 1.0, 0.2, 8.0, 260.0, 25.0] # upper bound

# optimisation
(xopt,fopt)  = optimise_MLE(fun,θG,lb,ub)




################### Plots for figures #################

# storing data for plots
t1_smooth = LinRange(0,11,1001)
return_smooth = twopoptwocompartmentreducedGreenspan(t1_smooth,xopt);

t1_smooth_plot = return_smooth[:,1]
ymle_smooth_tmp  = return_smooth[:,2]

len_t1_smooth_plot = length(t1_smooth_plot);
R_smooth_plot = zeros(len_t1_smooth_plot,1);
R1_smooth_plot = zeros(len_t1_smooth_plot,1);
Vp1c1_smooth_plot = zeros(len_t1_smooth_plot,1);
Vp2c1_smooth_plot = zeros(len_t1_smooth_plot,1);
Vp1c2_smooth_plot = zeros(len_t1_smooth_plot,1);
Vp2c2_smooth_plot = zeros(len_t1_smooth_plot,1);
Vp1c3_smooth_plot = zeros(len_t1_smooth_plot,1);
Vp2c3_smooth_plot = zeros(len_t1_smooth_plot,1);

for i=1:len_t1_smooth_plot
    R_smooth_plot[i] = ymle_smooth_tmp[i][9];
    R1_smooth_plot[i] = ymle_smooth_tmp[i][7];
    Vp1c1_smooth_plot[i] = ymle_smooth_tmp[i][1];
    Vp2c1_smooth_plot[i] = ymle_smooth_tmp[i][2];
    Vp1c2_smooth_plot[i] = ymle_smooth_tmp[i][3];
    Vp2c2_smooth_plot[i] = ymle_smooth_tmp[i][4];
    Vp1c3_smooth_plot[i] = ymle_smooth_tmp[i][5];
    Vp2c3_smooth_plot[i] = ymle_smooth_tmp[i][6];
end

index_to_lookup_endphasei = findfirst(R1_smooth_plot .> 0)[1] - 1;

t_end_phaseii = t1_smooth_plot[index_to_lookup_endphasei];


# plots 

day_shift = 2;
day_shift_t1_smooth_plot = day_shift*ones(len_t1_smooth_plot,1);

R_col = :green3;
R1_col = :black;

# Temporal evolution of structure
fig_pRR1R2_tmp = plot(plot_data,day_shift_t1_smooth_plot  + t1_smooth_plot,R_smooth_plot,c=R_col,lw=4,xlab=L"t \mathrm{ \ [days]}",ylab=L"R(t) \ [\mu \mathrm{m}]",legend=false,titlefont=fnt, guidefont=fnt, tickfont=fnt) # R(t)
fig_pRR1R2_tmp = plot!(fig_pRR1R2_tmp,day_shift_t1_smooth_plot  + t1_smooth_plot,R1_smooth_plot,c=R1_col,lw=4) # R1(t)
fig_pRR1R2_tmp = vline!(fig_pRR1R2_tmp,[(day_shift + t_end_phaseii)],c=:black,lw=2,linestyle=:dash)
display(fig_pRR1R2_tmp)
savefig(fig_pRR1R2_tmp,filepath_save[1] * "fig_pRR1R2.pdf")


# Temporal evolution of compartment 1, 2, 3 composition
fig_fracs_tmp = plot(day_shift_t1_smooth_plot  + t1_smooth_plot,Vp1c1_smooth_plot./(Vp1c1_smooth_plot + Vp2c1_smooth_plot),c=R_col,lw=4,xlim=(0,11),ylim=(0,1.1),xticks=[0,2,4,6,8,10],xlab=L"t \mathrm{ \ [days]}",ylab=L"\phi_{m}^{(j)}(t)\ [-]",legend=false,titlefont=fnt, guidefont=fnt, tickfont=fnt) # R(t)
fig_fracs_tmp = plot!(fig_fracs_tmp,day_shift_t1_smooth_plot  + t1_smooth_plot,Vp1c2_smooth_plot./(Vp1c2_smooth_plot + Vp2c2_smooth_plot),c=R1_col,lw=4,linestyle=:dash) # R1(t)
fig_fracs_tmp = vline!(fig_fracs_tmp,[(day_shift + t_end_phaseii)],c=:black,lw=2,linestyle=:dash)
display(fig_fracs_tmp)
savefig(fig_fracs_tmp,filepath_save[1] * "fig_fracs.pdf")

