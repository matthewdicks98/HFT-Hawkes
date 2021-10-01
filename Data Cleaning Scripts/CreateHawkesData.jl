using CSV, DataFrames, Dates, ProgressMeter, Plots, LaTeXStrings, TimeSeries, GR, JLD

# set working directory and clear the console
cd("C:\\Users\\Matt\\Desktop\\UCT Advanced Analytics\\HFT\\HFT-Lab-Exam")
clearconsole()

# function to create the Hawkes event type and time data
# output is a jld database
# only for 4 event types [MOB, MOS, LOB, LOS]
# assumes trades are from market orders, so no matching of limit orders
# this is also only for 1 day
function createHawkesData(date::String, ticker::String, write::Bool)

    full_data = CSV.read("test_data\\Clean\\TAQ\\JSECLEANTAQ"*ticker*"_c3.csv", DataFrame)
    println("Read in full data...")

    # get data only for a specific date
    formated_date = Date(date, "dd/mm/yyyy")
    data = full_data[findall(x -> x == formated_date, full_data[:,:date]), :]

    # want to create a dictionary where the keys are the events and the values are a list
    # of trade times where that event occured
    events_dict = Dict{String, Vector{Int64}}()

    # get all the trades
    trades = data[findall(x -> x == "TRADE", data[:,:eventType]), :]

    # get all market order buys
    mob = trades[findall(x -> x == 1, trades[:,:tradeSign]), :]

    # get all the market order sells
    mos = trades[findall(x -> x == -1, trades[:,:tradeSign]), :]

    # get all the limit order buy (bids)
    lob = data[findall(x -> x == "BID", data[:,:eventType]), :]

    # get all the limit order sells (asks)
    los = data[findall(x -> x == "ASK", data[:,:eventType]), :]

    # push event times to dict
    # convert times to all be in seconds, convert from nanoseconds to seconds

    # time zero is 9:00:00 need to reduce all other times to reflect amount of seconds from 0
    t_0 = Int64(9 * 3600)
    push!(events_dict, "MOB" => Int64.(Dates.value.(mob[:,:time])./1000000000 .- t_0))
    push!(events_dict, "MOS" => Int64.(Dates.value.(mos[:,:time])./1000000000 .- t_0))
    push!(events_dict, "LOB" => Int64.(Dates.value.(lob[:,:time])./1000000000 .- t_0))
    push!(events_dict, "LOS" => Int64.(Dates.value.(los[:,:time])./1000000000 .- t_0))

    # write to jld databse
    if write

        save("test_data\\Clean\\HAWKES\\HAWKES_"*date[1:2]*"_"*date[4:5]*"_"*date[7:end]*"_"*ticker*".jld", "HawkesData", events_dict)

    end

end

createHawkesData("10/07/2019","NPN", false)
createHawkesData("08/07/2019","NPN", false)
createHawkesData("08/07/2019","AGL", false)
