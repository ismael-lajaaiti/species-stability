using CSV
using DataFrames

df = DataFrame(CSV.File("data/white2020.csv"))

# Stack algae columns to ease data manipulation.
algae_names = names(df)[8:end]
df = stack(df, algae_names; variable_name = :species, value_name = :cover)

# Clean column names.
remove_trailing_whitespace(string) = replace(string, r"\s+$" => "")
new_colnames = remove_trailing_whitespace.(names(df))
new_colnames = lowercase.(new_colnames)
rename!(df, new_colnames)

# Missing value (only one) is set to 0.
df[ismissing.(df.cover), :cover] .= 0

# Convert the disturbance column to boolean values.
df.disturbance = [x == "Yes" for x in df.disturbance]

# Remove the treatment column as it contains redundant information.
# as well as the month column.
select!(df, Not(:treatment, :month))

# We also improving the naming of column for sake of clarity.
rename!(df, :diversity => :treatment)
rename!(df, :time => :month)

df.treatment = replace(
    remove_trailing_whitespace.(df.treatment),
    "+PGL" => :pgl,
    "-PGL" => :none,
    "-L" => :pg,
    "-G" => :pl,
    "-P" => :gl,
    "CONTROL" => :control,
)
df.treatment = Symbol.(df.treatment)

CSV.write("data/white2020_processed.csv", df)
