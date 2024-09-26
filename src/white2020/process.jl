using CSV
using DataFrames

df = DataFrame(CSV.File("data/white2020.csv"))
remove_trailing_whitespace(string) = replace(string, r"\s+$" => "")
new_colnames = remove_trailing_whitespace.(names(df))
rename!(df, new_colnames)

CSV.write("data/white2020_processed.csv", df)
