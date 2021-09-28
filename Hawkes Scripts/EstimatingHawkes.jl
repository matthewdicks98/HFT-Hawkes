# used to estimate a 4 dimensional hawkes process

using CSV, DataFrames, Dates, ProgressMeter, Plots, LaTeXStrings, TimeSeries, GR, JLD, Optim, StatsBase, StatsPlots, Distributions, HypothesisTests, LinearAlgebra

# set working directory and clear the console
cd("C:\\Users\\Matt\\Desktop\\UCT Advanced Analytics\\HFT\\HFT-Lab-Exam")
clearconsole()

# create Recursive function to be used in likelihood, computed using dynamic programming (more efficient)
function R(m, n, beta_m, events_times)

    # get the different events
    events = collect(keys(events_times))

    # set the beta
    beta_mn = beta_m[n]

    # get the event times for events m and n (for event m need first possible event to be 0)
    event_times_m = append!([0.0],events_times[events[m]])
    event_times_n = events_times[events[n]]

    # number of time steps that needs to be computed
    L = length(event_times_m)

    # initialize R_mn, R_mn(0) = 0 => R_mn(1) = 0 in code
    R_mn = zeros(Real, L+1, 1)

    # store the time step index of nth event that was just out of range for m != n
    prev_step = 1

    # loop over all the L times in the mth event
    for l in 2:L

        if n == m # self-exciting

            R_mn[l] = exp(-beta_mn * (event_times_m[l] - event_times_m[l-1])) * (1 + R_mn[l-1])

        else # mutually exciting

            left_term = exp(-beta_mn * (event_times_m[l] - event_times_m[l-1])) * R_mn[l-1]
            right_term = 0

            for k in prev_step:length(event_times_n)

                if event_times_n[k] >= event_times_m[l-1]

                    if event_times_n[k] < event_times_m[l]

                        right_term = right_term + exp(-beta_mn * (event_times_m[l] - event_times_n[k]))

                    else

                        prev_step = k
                        break # break the for loop if there all the other times are after the lth time of the mth event

                    end

                end

            end

            R_mn[l] = left_term + right_term

        end

    end

    return R_mn[2:L] # return all the recursive terms for all L time steps for nth event acting on mth

end

#R(1, 1, fill(0.2,4), events_times)[2]

# create the integrated intensity function for the mth event
function IntergratedIntensity_M(T, m, lambda_0_m, alpha_m, beta_m, events_times)

    # get the different events
    events = collect(keys(events_times))

    # base rate thoughout the window
    Lambda_m = lambda_0_m * T

    # loop over the events
    for n in 1:length(events)

        # get the exiting and decaying params for the nth event effecting the mth event
        alpha_mn = alpha_m[n]
        beta_mn = beta_m[n]

        # set the mth and nth event times
        event_times_m = events_times[events[m]]
        event_times_n = events_times[events[n]]

        # loop over the times in the nth event
        for i in 1:length(event_times_n)

            if event_times_n[i] <= T

                Lambda_m = Lambda_m + (alpha_mn/beta_mn) * (1 - exp(-beta_mn * (T - event_times_n[i])))

            end

        end

    end

    return Lambda_m

end

# create partial log likelihood
function PartialLogLikelihood(T, m, lambda_0_m, alpha_m, beta_m, events_times)

    # get the different events
    events = collect(keys(events_times))

    # get the mth event data
    event_times_m = events_times[events[m]]

    # set the first 2 terms of the sum
    ll_m = T - IntergratedIntensity_M(T, m, lambda_0_m, alpha_m, beta_m, events_times)

    # precompute R matrix
    R_m = zeros(Real, length(events), length(event_times_m))
    for n in 1:length(events)

        # get all the recursive values for event n acting on event m
        R_m[n,:] = R(m, n, beta_m, events_times)

    end

    # loop over all the mth event times
    for l in 1:length(event_times_m)

        temp = lambda_0_m

        for n in 1:length(events)

            temp = temp + alpha_m[n] * R_m[n,l]

        end

        if temp > 0 # just make sure positive

            ll_m = ll_m + log(temp)

        end

    end

    return ll_m

end

#PartialLogLikelihood(60589 - 32450, 3, 0.01, fill(0.01,4), fill(0.2,4), events_times)
#5962

# create full log likelihood (this is the function to be optimized)
function LogLikelihood(params, T, events_times)

    # get the different events
    events = collect(keys(events_times))

    # get number of events
    M = length(collect(keys(events_times)))

    # get params in the correct format to input into the likelihood
    lambda_0 = params[1:M]
    alpha = reshape(params[(M+1):((M * M) + M)], M, M)
    beta = reshape(params[((M * M) + M + 1):end], M,  M)

    # initialize the log likelihood
    ll = 0

    # loop over all the event types, compute the partial and sum them to get the full
    for m in 1:length(events)

        ll = ll + PartialLogLikelihood(T, m, lambda_0[m], alpha[m,:], beta[m,:], events_times)

    end

    return -ll

end

function IntergratedDensity_M(T2, T1, m, lambda_0_m, alpha_m, beta_m, events_times)

    # get the different events
    events = collect(keys(events_times))

    event_times_m = events_times[events[m]]

    # precompute R
    R_m = zeros(Real, length(events), length(event_times_m))
    for n in 1:length(events)

        # get all the recursive values for event n acting on event m
        R_m[n,:] = R(m, n, beta_m, events_times)

    end

    # base rate thoughout the window
    Lambda_m = lambda_0_m * (T2 - T1)

    # loop over the events
    for n in 1:length(events)

        # get the exiting and decaying params for the nth event effecting the mth event
        alpha_mn = alpha_m[n]
        beta_mn = beta_m[n]

        # set the mth and nth event times
        event_times_m = events_times[events[m]]
        event_times_n = events_times[events[n]]

        Lambda_m = Lambda_m + (alpha_mn/beta_mn) * R_m[n,findall(x -> x == T1, event_times_m)[1]] * (1 - exp(-beta_mn * (T2 - T1)))

        for k in 1:length(event_times_n)

            if event_times_n[k] >= T1 && event_times_n[k] < T2

                Lambda_m = Lambda_m + (1 - exp(-beta_mn * (T2 - event_times_n[k])))

            end

        end

    end

    return Lambda_m

end

function Durations(params_hat, events_times)

    # get the different events
    events = collect(keys(events_times))

    # get number of events
    M = length(collect(keys(events_times)))

    # get params in the correct format to input into the likelihood
    lambda_0 = params_hat[1:M]
    alpha = reshape(params_hat[(M+1):((M * M) + M)], M, M)
    beta = reshape(params_hat[((M * M) + M + 1):end], M,  M)

    durations = []

    println("Starting Compute")

    @showprogress "Computing Durations..." for m in 1:M

        # get the times for the mth event
        event_times_m = events_times[events[m]]

        # create storage for the durations
        intensities = fill(0.0, length(event_times_m))

        # get the mth event params
        lambda_0_m = lambda_0[m]
        alpha_m = alpha[m,:]
        beta_m = beta[m,:]

        for i in 1:length(event_times_m)

            intensities[i] = IntergratedIntensity_M(event_times_m[i], m, lambda_0_m, alpha_m, beta_m, events_times)

        end

        push!(durations, diff(intensities))

        #for i in 1:length(event_times_m)

        #    intensities[i] = IntergratedDensity_M(event_times_m[i], event_times_m[i-1], m, lambda_0_m, alpha_m, beta_m, events_times)

        #end

        #push!(durations, intensities)

    end

    return durations

end

function SpectralRadius(params_hat, events_times)

    # get the different events
    events = collect(keys(events_times))

    # get number of events
    M = length(collect(keys(events_times)))

    # get params in the correct format to input into the likelihood
    lambda_0 = params_hat[1:M]
    alpha = reshape(params_hat[(M+1):((M * M) + M)], M, M)
    beta = reshape(params_hat[((M * M) + M + 1):end], M,  M)

    Gamma = fill(0.0, M, M)

    for m in 1:M

        for n in 1:M

            Gamma[m,n] = alpha[m,n] / beta[m,n]

        end

    end

    sr = maximum(abs.(eigvals(Gamma)))

    return sr

end

function PlotDurations(durations, events_times, write)

    # get the different events
    events = collect(keys(events_times))

    # plot first
    q1 = qqplot(Exponential(1), validation_durations[1], color = 1, label = collect(keys(events_times))[1], legend = :outertopright,
    markercolor = 1, markerstrokecolor = 1, markershape=:circle, title = "Durations QQ-plot for "*events[1]*" - Exponential(1)",
    xlabel = "Theoretical Quantiles", ylabel = "Sample Quantiles")
    q2 = qqplot(Exponential(1), validation_durations[2], color = 2, label = collect(keys(events_times))[2], legend = :outertopright,
    markercolor = 2, markerstrokecolor = 2, markershape=:circle, title = "Durations QQ-plot for "*events[2]*" - Exponential(1)",
    xlabel = "Theoretical Quantiles", ylabel = "Sample Quantiles")
    q3 = qqplot(Exponential(1), validation_durations[3], color = 3, label = collect(keys(events_times))[3], legend = :outertopright,
    markercolor = 3, markerstrokecolor = 3, markershape=:circle, title = "Durations QQ-plot for "*events[3]*" - Exponential(1)",
    xlabel = "Theoretical Quantiles", ylabel = "Sample Quantiles")
    q4 = qqplot(Exponential(1), validation_durations[4], color = 4, label = collect(keys(events_times))[4], legend = :outertopright,
    markercolor = 4, markerstrokecolor = 4, markershape=:circle, title = "Durations QQ-plot for "*events[4]*" - Exponential(1)",
    xlabel = "Theoretical Quantiles", ylabel = "Sample Quantiles")

    if write

        Plots.savefig(q1, "Images\\QQ_"*events[1]*"_NPN.pdf")
        Plots.savefig(q2, "Images\\QQ_"*events[2]*"_NPN.pdf")
        Plots.savefig(q3, "Images\\QQ_"*events[3]*"_NPN.pdf")
        Plots.savefig(q4, "Images\\QQ_"*events[4]*"_NPN.pdf")

    end

    display(q1)
    display(q2)
    display(q3)
    display(q4)

end

function EventVisualization(params_hat, m, events_times, write)

    # get the different events
    events = collect(keys(events_times))

    # get number of events
    M = length(collect(keys(events_times)))

    # get params in the correct format to input into the likelihood
    lambda_0 = params_hat[1:M]
    alpha = reshape(params_hat[(M+1):((M * M) + M)], M, M)
    beta = reshape(params_hat[((M * M) + M + 1):end], M,  M)

    # get branching ratios for event type m
    Gamma_mn = alpha[:,m] ./ beta[:,m]

    # get the half lives of the excitation event alpha_mn
    half_lives = log(2) ./ beta[:,m]

    market_scale = max(length(events_times[events[2]]), length(events_times[events[4]]))/10
    limit_scale = max(length(events_times[events[1]]), length(events_times[events[3]]))/20

    s = Plots.scatter([half_lives[1]], [Gamma_mn[1]], color = :red, label = events[1], legend = :outertopright, markercolor = 1,
    markerstrokecolor = 1, markershape=:circle, markersize = Int(round(length(events_times[events[1]])/limit_scale)),
    title = "Mutual and Self-excitations due to "*events[1]*"", xlabel = "Half-life", ylabel = "# Excitations", ylim = (-0.1, 1))

    for n in 2:length(events)

        if n % 2 == 0

            Plots.scatter!([half_lives[n]], [Gamma_mn[n]], color = :red, label = events[n], legend = :outertopright, markercolor = n,
             markerstrokecolor = n, markershape=:circle, markersize = Int(round(length(events_times[events[n]])/market_scale)),
             title = "Mutual and Self-excitations due to "*events[n]*"", xlabel = "Half-life", ylabel = "# Excitations", ylim = (-0.1, 1))

        else

            Plots.scatter!([half_lives[n]], [Gamma_mn[n]], color = :red, label = events[n], legend = :outertopright, markercolor = n,
            markerstrokecolor = n, markershape=:circle, markersize = round(length(events_times[events[n]])/limit_scale),
            title = "Mutual and Self-excitations due to "*events[n]*"", xlabel = "Half-life", ylabel = "# Excitations", ylim = (-0.1, 1))

        end

    end

    if write

        Plots.savefig(s, "Images\\Excitations_"*events[m]*"_NPN.pdf")

    end

    display(s)

end

function PlotACF(durations, lags, events_times, write)

    # get the different events
    events = collect(keys(events_times))

    a1 = Plots.plot(1:lags, autocor(durations[1], 1:lags), seriestype = :sticks, color = :black, title = "Auto-correlation in Durations - "*events[1]*"",
    xlabel = "Lags", ylabel = "Auto-correlation")
    Plots.hline!(1:lags, [1.96/sqrt(length(durations[1]))], seriestype = :line, color = :red, legend = false, lw = 2)
    Plots.hline!(1:lags, [-1.96/sqrt(length(durations[1]))], seriestype = :line, color = :red, legend = false, lw = 2)

    a2 = Plots.plot(1:lags, autocor(durations[2], 1:lags), seriestype = :sticks, color = :black, title = "Auto-correlation in Durations - "*events[2]*"",
    xlabel = "Lags", ylabel = "Auto-correlation")
    Plots.hline!(1:lags, [1.96/sqrt(length(durations[2]))], seriestype = :line, color = :red, legend = false, lw = 2)
    Plots.hline!(1:lags, [-1.96/sqrt(length(durations[2]))], seriestype = :line, color = :red, legend = false, lw = 2)

    a3 = Plots.plot(1:lags, autocor(durations[3], 1:lags), seriestype = :sticks, color = :black, title = "Auto-correlation in Durations - "*events[3]*"",
    xlabel = "Lags", ylabel = "Auto-correlation")
    Plots.hline!(1:lags, [1.96/sqrt(length(durations[3]))], seriestype = :line, color = :red, legend = false, lw = 2)
    Plots.hline!(1:lags, [-1.96/sqrt(length(durations[3]))], seriestype = :line, color = :red, legend = false, lw = 2)

    a4 = Plots.plot(1:lags, autocor(durations[4], 1:lags), seriestype = :sticks, color = :black, title = "Auto-correlation in Durations - "*events[4]*"",
    xlabel = "Lags", ylabel = "Auto-correlation")
    Plots.hline!(1:lags, [1.96/sqrt(length(durations[4]))], seriestype = :line, color = :red, legend = false, lw = 2)
    Plots.hline!(1:lags, [-1.96/sqrt(length(durations[4]))], seriestype = :line, color = :red, legend = false, lw = 2)

    if write

        Plots.savefig(a1, "Images\\ACF_"*events[1]*"_NPN.pdf")
        Plots.savefig(a2, "Images\\ACF_"*events[2]*"_NPN.pdf")
        Plots.savefig(a3, "Images\\ACF_"*events[3]*"_NPN.pdf")
        Plots.savefig(a4, "Images\\ACF_"*events[4]*"_NPN.pdf")

    end

    display(a1)
    display(a2)
    display(a3)
    display(a4)

end

# get the hawkes data, the times each event type occurs
events_times = load("test_data\\Clean\\HAWKES\\HAWKES_10_07_2019_NPN.jld")["HawkesData"]

# set starting params
events = collect(keys(events_times))
M = length(collect(keys(events_times)))
lambda_0_initial = [0.2, 0.2, 0.2, 0.2]
alpha_initial = fill(0.01, M*M)
beta_initial = fill(0.2, M*M)
T = Int64((16 * 3600 + 49 * 60 + 49) - (9 * 3600))

params_initial = append!(append!(lambda_0_initial, alpha_initial), beta_initial)

# perform the optimization
lower = fill(0.0, length(params_initial)) # params >= 0
upper = fill(Inf, length(params_initial))

#fun = TwiceDifferentiable(params -> LogLikelihood(params, T, events_times), params_initial)
#fun_constraints = TwiceDifferentiableConstraints(lower, upper)

#@time res1 = optimize(params -> LogLikelihood(params, T, events_times), params_initial, LBFGS(), autodiff = :forward, Optim.Options(show_trace = true, iterations = 10))

@time res = optimize(params -> LogLikelihood(params, T, events_times), lower, upper, params_initial, Fminbox(LBFGS()), autodiff = :forward, Optim.Options(store_trace = true, show_trace = true, outer_iterations = 5, iterations = 200))

println(res)
println(Optim.converged(res))
println(Optim.minimizer(res))
Plots.plot(1:length(Optim.f_trace(res)), Optim.f_trace(res), seriestype = :line)

# compute and plot the durations wrt the Exp(1) distribution
validation_durations = Durations(Optim.minimizer(res), events_times)

# Exponential check
PlotDurations(validation_durations, events_times, false)

# acf check for stationarity
PlotACF(validation_durations, 100, events_times, false)

# Spectral radius check for stationarity
println(SpectralRadius(Optim.minimizer(res), events_times))

# create the event generations Plots
EventVisualization(Optim.minimizer(res), 1, events_times, false)
EventVisualization(Optim.minimizer(res), 2, events_times, false)
EventVisualization(Optim.minimizer(res), 3, events_times, false)
EventVisualization(Optim.minimizer(res), 4, events_times, false)

# Hypothesis tests to confirm iid exp durations

# Exp assumptions
#println(round(pvalue(ExactOneSampleKSTest(validation_durations[1], Exponential(1))), digits = 4))
#println(round(pvalue(ExactOneSampleKSTest(validation_durations[2], Exponential(1))), digits = 4))
#println(round(pvalue(ExactOneSampleKSTest(validation_durations[3], Exponential(1))), digits = 4))
#println(round(pvalue(ExactOneSampleKSTest(validation_durations[4], Exponential(1))), digits = 4))

# iid assumtions
#println(round(pvalue(LjungBoxTest(validation_durations[1], 100, 1)), digits = 4))
#println(round(pvalue(LjungBoxTest(validation_durations[2], 100, 1)), digits = 4))
#println(round(pvalue(LjungBoxTest(validation_durations[3], 100, 1)), digits = 4))
#println(round(pvalue(LjungBoxTest(validation_durations[4], 100, 1)), digits = 4))
