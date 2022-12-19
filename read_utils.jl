function read_price_data(filename)
    half_hour_df = DataFrame(XLSX.readtable(filename, "Half-hourly data"))
    daily_df = DataFrame(XLSX.readtable(filename, "Daily data"))

    market_3_prices = zeros(nrow(half_hour_df))
    for d in 1:nrow(daily_df)
        market_3_prices[(d - 1) * 48 + 1:d * 48] .= daily_df[d, "Market 3 Price [£/MWh]"]
    end

    price_data = DataFrame(
        time = DateTime.(half_hour_df[:, "Time"]),
        market_1 = half_hour_df[:, "Market 1 Price [£/MWh]"],
        market_2 = half_hour_df[:, "Market 2 Price [£/MWh]"],
        market_3 = market_3_prices
        )

    return price_data
end

function read_battery_parameters(filename)
    battery_parameters_df = DataFrame(XLSX.readtable(filename, "Data"))
    battery_parameters = Dict(
        battery_parameters_df[i, "Parameters"] => battery_parameters_df[i, "Values"]
        for i in 1:nrow(battery_parameters_df)
        )

    return battery_parameters
end
