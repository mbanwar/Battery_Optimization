function optimize_battery_operation(
    battery_parameters,
    all_price_data,
    cycling_on,
    horizon_start,
    simulation_months)

    # Define markets and granularity
    markets = ["market_1", "market_2", "market_3"]

    # Initialize storage level
    init_storage = 0.0
    init_discharge = 0.0

    # Battery parameters parsing
    maximum_discharge_rate = battery_parameters["Max discharging rate"]
    maximum_charge_rate = battery_parameters["Max charging rate"]
    charging_efficiency = 1 - battery_parameters["Battery charging efficiency"]
    discharging_efficiency = 1 - battery_parameters["Battery discharging efficiency"]
    max_storage_level = battery_parameters["Max storage volume"]
    min_storage_level = 0.0 # Can be set to any non-zero value to avoid deep discharge
    max_cycles = battery_parameters["Lifetime (2)"]
    battery_degradation = battery_parameters["Storage volume degradation rate"] * 1e-2
    fixed_OnM = battery_parameters["Fixed Operational Costs"]

    "Output File"
    output_dir = joinpath(@__DIR__, "results")
    mkpath(output_dir)
    output_file = joinpath(output_dir, "results_cycling_$(cycling_on).xlsx")

    # Output Data
    p_in_data = DataFrame(time = [], market_1 =[], market_2 = [], market_3 = [])
    p_out_data = DataFrame(time = [], market_1 =[], market_2 = [], market_3 = [])
    profit_data = DataFrame(time = [], market_1 =[], market_2 = [], market_3 = [])
    storage_level_data = DataFrame(time = [], storage = [])
    cycle_data = DataFrame(date = [], number_cycles = [])

    delta_t = 0.5 # 30-min time interval

    for m in 1:simulation_months
        start_date = horizon_start + Month(m - 1)
        end_date = horizon_start + Month(m) - Day(1)
        println("Start Date: ", start_date)
        println("End Date: ", end_date)

        timestamps = collect(DateTime(start_date):Minute(30):DateTime(end_date) + Hour(23) + Minute(30))
        time_horizon = collect(1:length(timestamps))

        price_data = filter(row -> row.time in timestamps, all_price_data)
        # JuMP Model
        model = Model(GLPK.Optimizer)

        # Add decision variables
        @variable(model, p_in[m in markets, t in time_horizon] >= 0)
        @variable(model, p_out[m in markets, t in time_horizon] >= 0)

        # Variables for modelling battery cycling
        @variable(model, storage_inc[t in time_horizon] >= 0)
        @variable(model, storage_dec[t in time_horizon] >= 0)
        @variable(model, dummy == 0.0)

        # Add constraints and expressions
        @info "Adding Constraints"

        # Storage evolution expression
        @expression(model, storage_level[t in time_horizon],
            delta_t * sum(p_in[m, t] * charging_efficiency - p_out[m, t] / discharging_efficiency for m in markets))

        # Expressions for calculating storage decrease at each timestep.
        @expression(model, delta_storage[t in time_horizon], dummy * 1.0)
        @expression(model, discharge_sum[t in time_horizon], storage_dec[t] * 1.0)
        @expression(model, num_cycles[t in time_horizon], dummy * 1.0)

        for t in time_horizon
            if t == first(time_horizon)
                add_to_expression!(storage_level[t], init_storage)
                add_to_expression!(delta_storage[t], storage_level[t] - init_storage)
                add_to_expression!(discharge_sum[t], init_discharge)
            else
                add_to_expression!(storage_level[t], storage_level[t - 1])
                add_to_expression!(delta_storage[t], storage_level[t] - storage_level[t - 1])
                add_to_expression!(discharge_sum[t], discharge_sum[t - 1])

                # Constant daily trading for market 3
                if Date(timestamps[t]) == Date(timestamps[t - 1])
                    @constraint(model, p_in["market_3", t] == p_in["market_3", t - 1])
                    @constraint(model, p_out["market_3", t] == p_out["market_3", t - 1])
                end
            end

            # Assuming 1 discharge cycle corresponds to once the total discharge sum reaches the maximum storage level
            add_to_expression!(num_cycles[t], discharge_sum[t] / max_storage_level)

            # Maximum charging rate constraint
            @constraint(model, sum(p_in[m, t] for m in markets) <= maximum_charge_rate)

            # Maximum discharging rate constraint
            @constraint(model, sum(p_out[m, t] for m in markets) <= maximum_discharge_rate)

            # Storage Level constraints
            @constraint(model, storage_level[t] >= min_storage_level)
            @constraint(model, storage_level[t] <= max_storage_level * (1 - battery_degradation * num_cycles[t] * cycling_on))

            # Storage increase and decrease modelling
            @constraint(model, delta_storage[t] == storage_inc[t] - storage_dec[t])

            # Limiting num_cycles
            @constraint(model, num_cycles[t] <= max_cycles)
        end


        cycle_scaling = 1e-3 # To weigh the storage inc and dec accuracy
        @objective(model, Max,
            sum(sum((p_out[m, t] - p_in[m, t]) * price_data[t, m] for m in markets) -
                storage_inc[t]* cycling_on * cycle_scaling -
                storage_dec[t]* cycling_on * cycle_scaling for t in time_horizon) * delta_t)


        @info "Solving Model"
        @time optimize!(model)
        println("Solver Status: ", termination_status(model))
        println("Objective Value: ", objective_value(model))

        # Store output data
        for t in time_horizon
            push!(p_in_data, [timestamps[t] value(p_in["market_1", t]) value(p_in["market_2", t]) value(p_in["market_3", t])])
            push!(p_out_data, [timestamps[t] value(p_out["market_1", t]) value(p_out["market_2", t]) value(p_out["market_3", t])])
            push!(storage_level_data, [timestamps[t] value(storage_level[t])])
            interval_profit = [(value(p_out[m, t]) - value(p_in[m, t])) * delta_t * price_data[t, m] for m in markets]'
            push!(profit_data, [timestamps[t] interval_profit])
        end

        # Update storage level data for rolling optimization
        init_storage = last(storage_level_data[:, "storage"])
        init_discharge = last(value.(discharge_sum))

        if cycling_on
            push!(cycle_data, [end_date last(value.(num_cycles))])
        end
    end

    # Calculate yearly profit
    yearly_profit_data = DataFrame(year = [], profit = [])
    for y in unique(Year.(profit_data[:, "time"]))
        profit = sum(profit_data[t, m] for m in markets, t in 1:nrow(profit_data) if Year(profit_data[t, "time"]) == y)
        profit -= fixed_OnM
        push!(yearly_profit_data, [Dates.value(y) profit])
    end


    XLSX.writetable(output_file,
    "Charging" => p_in_data,
    "Discharging" => p_out_data,
    "Storage_Level" => storage_level_data,
    "Profit_Timeseries" => profit_data,
    "Yearly_Profit" => yearly_profit_data,
    "Number_of_Cycles" => cycle_data,
    overwrite = true)

end
