using DataFrames
using Dates
using GLPK
using JuMP
using XLSX

include("battery_optimization.jl")
include("read_utils.jl")

# Battery specifications
battery_parameters = read_battery_parameters(joinpath("data", "Battery Parameters.xlsx"))

# Whether battery cycling is modelled
cycling_on = false

# Read Price Data
all_price_data = read_price_data(joinpath("data", "Market Data.xlsx"))

# Time indices
horizon_start = Date(2018, 01, 01)
simulation_months = 36

optimize_battery_operation(
    battery_parameters,
    all_price_data,
    cycling_on,
    horizon_start,
    simulation_months)
