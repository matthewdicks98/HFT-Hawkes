# simulation of a 4 dimensional Hawkes process
using CSV, DataFrames, Dates, ProgressMeter, Plots, LaTeXStrings, TimeSeries, GR, JLD, Optim, StatsBase, StatsPlots, Distributions, HypothesisTests, LinearAlgebra, Random

# set working directory and clear the console
cd("C:\\Users\\Matt\\Desktop\\UCT Advanced Analytics\\HFT\\HFT-Lab-Exam")
clearconsole()

function Lambda(m, t, lambda_0, alpha, beta, events_times)

    # get the different events
    events = collect(keys(events_times))

    # fetch the params
    lambda_0_m = lambda_0[m]
    alpha_m = alpha[m,:]
    beta_m = beta[m,:]

    # set the initial value for lambda which is at least lambda_0_m
    Lambda_m = lambda_0_m

    # compute the increase in intensity from the other events
    for n in 1:length(events)

        # event times for the event n
        event_times_n = events_times[events[n]]

        # set the correct params
        alpha_mn = alpha_m[n]
        beta_mn = beta_m[n]

        for i in 1:length(event_times_n)

            if event_times_n[i] < t

                Lambda_m = Lambda_m + alpha_mn * exp(-beta_mn * (t - event_times_n[i]))

            else

                break

            end

        end

    end

    return Lambda_m

end

function Compute_I_M(t, lambda_0, alpha, beta, events_times)

    # get the different events
    events = collect(keys(events_times))

    # set storage
    intensities = fill(0.0, length(events))

    # add up all the internsities
    for m in 1:length(events)

        intensities[m] = Lambda(m, t, lambda_0, alpha, beta, events_times)

    end

    return intensities

end

function Find_n0(intensities, I_star, D, events_times)

    # get the different events
    events = collect(keys(events_times))

    I = append!([0.0], cumsum(intensities))

    #println(I)
    #println(I_star)
    #println(D)

    for n in 2:length(I)

        if D > I[n-1]/I_star && D <= I[n]/I_star

            return n - 1

        end

    end

end
intensities = Compute_I_M(20000, [0.2, 0.1, 0.3, 0.15], fill(0.01, 4, 4), fill(0.2, 4, 4), events_times)
Find_n0(intensities, sum(intensities), rand(Uniform(0,1)), events_times)

function SimulateHawkesProcess(T, lambda_0, alpha, beta, seed)

    println("Starting Simulation")

    # initialize events dict and combined times
    events = ["LOB", "MOB", "LOS", "MOS"]
    events_times = Dict()
    for i in 1:length(events)
        push!(events_times, events[i] => [])
    end
    combined_times = []

    # initialize counters
    i = 0
    i_m_counters = fill(0, length(events))

    # initialize intensities and max intensity
    intensities = lambda_0
    I_star = sum(lambda_0)

    # keep track of lat event type to have occured
    n0 = 0

    # generate first event
    #Random.seed!(seed)
    U = rand(Uniform(0,1))
    s = -(1/I_star) * log(U)
    if s <= T

        # find event type n0 for the time
        D = rand(Uniform(0,1))
        n0 = Find_n0(intensities, I_star, D, events_times)

        # add time to that event type and update counter
        push!(events_times[events[n0]], s)
        i_m_counters[n0] += 1

        # add time to total times and update counter
        push!(combined_times, events_times[events[n0]][end])
        i += 1

    else

        return events_times

    end

    # general routine
    while true

        # update I_star and intensities
        if length(combined_times) == 1 # if only 1 event has occured then prev intensity is just lambda 0

            I_star = sum(intensities) + sum(alpha[:,n0])

        else # more than one event has occured

            # get the previous time an event occured
            prev_time = combined_times[end-1]

            # compute intensity due to previous time step and add effect of current n0 event
            intensities = Compute_I_M(combined_times[end-1], lambda_0, alpha, beta, events_times)
            I_star = sum(intensities) + sum(alpha[:,n0])

        end

        seed += 1
        #Random.seed!(seed)

        while true

            # generate new exponential
            U = rand(Uniform(0, 1))
            s = s - (1/I_star) * log(U)

            seed += 1
            #Random.seed!(seed)
            if s <= T

                # find event type to for the time
                D = rand(Uniform(0,1))

                # find intensities up to time s
                I_M = sum(Compute_I_M(s, lambda_0, alpha, beta, events_times))

                if D <= (I_M/I_star)

                    # find event type to occur
                    n0 = Find_n0(Compute_I_M(s, lambda_0, alpha, beta, events_times), I_star, D, events_times)

                    # add time to that event type and update counter
                    push!(events_times[events[n0]], s)
                    i_m_counters[n0] += 1

                    # add time to total times and update counter
                    push!(combined_times, events_times[events[n0]][end])
                    i += 1

                    break # break out of while loop

                else

                    I_star = sum(Compute_I_M(s, lambda_0, alpha, beta, events_times))

                end

            else

                return events_times

            end

        end

    end

    return events_times

end

function Durations(lambda_0, alpha, beta, events_times)

    # get the different events
    events = collect(keys(events_times))

    # get number of events
    M = length(collect(keys(events_times)))

    # get params in the correct format to input into the likelihood
    #lambda_0 = params_hat[1:M]
    #alpha = reshape(params_hat[(M+1):((M * M) + M)], M, M)
    #beta = reshape(params_hat[((M * M) + M + 1):end], M,  M)

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

        Plots.savefig(q1, "Images\\QQ_"*events[1]*"_SIM.pdf")
        Plots.savefig(q2, "Images\\QQ_"*events[2]*"_SIM.pdf")
        Plots.savefig(q3, "Images\\QQ_"*events[3]*"_SIM.pdf")
        Plots.savefig(q4, "Images\\QQ_"*events[4]*"_SIM.pdf")

    end

    display(q1)
    display(q2)
    display(q3)
    display(q4)

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

        Plots.savefig(a1, "Images\\ACF_"*events[1]*"_SIM.pdf")
        Plots.savefig(a2, "Images\\ACF_"*events[2]*"_SIM.pdf")
        Plots.savefig(a3, "Images\\ACF_"*events[3]*"_SIM.pdf")
        Plots.savefig(a4, "Images\\ACF_"*events[4]*"_SIM.pdf")

    end

    display(a1)
    display(a2)
    display(a3)
    display(a4)

end

function PlotEvents(events_times, write)

    # set number of events for the least frequent counting process
    num_events = 100

    # get the different events
    events = collect(keys(events_times))

    # get least frequent event type (MOS) last time and use this as the max time in the window
    T = events_times[events[1]][num_events]

    events_counter = Dict(events[1] => [[],[]], events[2]  => [[],[]], events[3]  => [[],[]], events[4]  => [[],[]])

    # build lines
    for m in 1:length(events)

        # get the event times for mth event
        event_times_m = events_times[events[m]]

        # loop over all the events and add the times and count the events
        for i in 1:length(event_times_m)

            if event_times_m[i] <= T

                push!(events_counter[events[m]][2], event_times_m[i])

                # if there are no events then add one
                if events_counter[events[m]][1] == []

                    push!(events_counter[events[m]][1], 1)

                else # increment counter

                    push!(events_counter[events[m]][1], events_counter[events[m]][1][end] + 1)

                end

            else

                break

            end

        end

    end

    p1 = Plots.plot(events_counter[events[1]][2], events_counter[events[1]][1], linetype=:steppre,
    label = events[1], legend = :outertopright, markercolor = 1, markerstrokecolor = 1, markershape=:circle, markersize = 3,
    title = "Visualizing the Hawkes Process's Counts", xlabel = "Time", ylabel = "# Counts")
    Plots.plot!(events_counter[events[2]][2], events_counter[events[2]][1], linetype=:steppost,
    label = events[2], legend = :outertopright, markercolor = 2, markerstrokecolor = 2, markershape=:circle, markersize = 3)
    Plots.plot!(events_counter[events[3]][2], events_counter[events[3]][1], linetype=:steppost,
    label = events[3], legend = :outertopright, markercolor = 3, markerstrokecolor = 3, markershape=:circle, markersize = 3)
    Plots.plot!(events_counter[events[4]][2], events_counter[events[4]][1], linetype=:steppost,
    label = events[4], legend = :outertopright, markercolor = 4, markerstrokecolor = 4, markershape=:circle, markersize = 3)

    p2 = Plots.plot(events_counter[events[1]][2], fill(1,length(events_counter[events[1]][2])), linetype=:steppost,
    label = events[1], legend = :outertopright, markercolor = 1, markerstrokecolor = 1, markershape=:circle,
    title = "Visualizing the Hawkes Process's Events", xlabel = "Time", ylabel = "Event Types")
    Plots.plot!(events_counter[events[2]][2], fill(2,length(events_counter[events[2]][2])), linetype=:steppost,
    label = events[2], legend = :outertopright, markercolor = 2, markerstrokecolor = 2, markershape=:circle)
    Plots.plot!(events_counter[events[3]][2], fill(3,length(events_counter[events[3]][2])), linetype=:steppost,
    label = events[3], legend = :outertopright, markercolor = 3, markerstrokecolor = 3, markershape=:circle)
    Plots.plot!(events_counter[events[4]][2], fill(4,length(events_counter[events[4]][2])), linetype=:steppost,
    label = events[4], legend = :outertopright, markercolor = 4, markerstrokecolor = 4, markershape=:circle)

    if write

        Plots.savefig(p1, "Images\\COUNTS_SIM.pdf")
        Plots.savefig(p2, "Images\\EVENTS_SIM.pdf")

    end

    display(p1)
    display(p2)

end

params_initial = Float64[]
open("Params/Estimated_Parameters.txt", "r") do file

    for line in readlines(file)

        # print the line
        push!(params_initial, parse(Float64, line))

    end

end

# get number of events
M = 4

# get params in the correct format to input into the likelihood
lambda_0_initial = params_initial[1:M]
alpha_initial = reshape(params_initial[(M+1):((M * M) + M)], M, M)
beta_initial = reshape(params_initial[((M * M) + M + 1):end], M,  M)

# set params
#lambda_0_initial = [0.1, 0.2, 0.4, 0.8]
#alpha_initial = fill(0.01, M, M)
#beta_initial = fill(0.2, M, M)
T = Int64((16 * 3600 + 49 * 60 + 49) - (9 * 3600))

@time events_times = SimulateHawkesProcess(T, lambda_0_initial, alpha_initial, beta_initial, 1)

# validate the simulation by checking durations
validation_durations = Durations(lambda_0_initial, alpha_initial, beta_initial, events_times)

# Exponential check
PlotDurations(validation_durations, events_times, false)

# acf check for stationarity
PlotACF(validation_durations, 100, events_times, false)

# visualize the Hawkes Process
PlotEvents(events_times, false)

# if happy write event times to a jld database
save("test_data\\Clean\\HAWKES\\HAWKES_SIMULATED.jld", "HawkesDataSim", events_times)

# count sim
length(load("test_data\\Clean\\HAWKES\\HAWKES_SIMULATED.jld")["HawkesDataSim"]["LOB"])
length(load("test_data\\Clean\\HAWKES\\HAWKES_SIMULATED.jld")["HawkesDataSim"]["MOB"])
length(load("test_data\\Clean\\HAWKES\\HAWKES_SIMULATED.jld")["HawkesDataSim"]["LOS"])
length(load("test_data\\Clean\\HAWKES\\HAWKES_SIMULATED.jld")["HawkesDataSim"]["MOS"])

# count obs
length(load("test_data\\Clean\\HAWKES\\HAWKES_10_07_2019_NPN.jld")["HawkesData"]["LOB"])
length(load("test_data\\Clean\\HAWKES\\HAWKES_10_07_2019_NPN.jld")["HawkesData"]["MOB"])
length(load("test_data\\Clean\\HAWKES\\HAWKES_10_07_2019_NPN.jld")["HawkesData"]["LOS"])
length(load("test_data\\Clean\\HAWKES\\HAWKES_10_07_2019_NPN.jld")["HawkesData"]["MOS"])
