# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

#=
import numpy as np
import pandas as pd
import pyarrow as pa
from scipy.interpolate import BSpline

def read_evaluate_and_save_splines(input_file, output_file, x_eval):
    """
    Read spline coefficients, evaluate splines, and save results to Arrow file
    
    Parameters:
    input_file (str): Path to input Arrow file containing coefficients
    output_file (str): Path where output Arrow file will be saved
    x_eval (array-like): x values at which to evaluate the splines
    
    Returns:
    str: Path to the saved Arrow file
    """
    # Read the Arrow file
    arrow_table = pa.ipc.RecordBatchFileReader(
        pa.memory_map(input_file, 'r')
    ).read_all()
    df = arrow_table.to_pandas()
    
    # Define spline parameters
    knots = np.array([6.0, 13.0, 20.0, 27.0, 34.0, 41.0, 48.0, 55.0])
    degree = 3
    
    # Initialize list to store results
    spline_values = []
    
    # Evaluate splines for each set of coefficients
    for i in range(df.shape[0]):
        c = np.array(df.coefficients[i])
        spl = BSpline(knots, c, degree)
        y = spl(x_eval)
        spline_values.append(y)
    
    # Create a DataFrame with results
    results_df = pd.DataFrame({
        'annotation': df.annotation,
        'spline_values': spline_values,
        'x_values': [x_eval] * len(df)  # Store x values for each spline
    })
    
    # Convert DataFrame to Arrow Table
    results_table = pa.Table.from_pandas(results_df)
    
    # Write to Arrow file
    with pa.OSFile(output_file, 'wb') as sink:
        writer = pa.ipc.new_file(sink, results_table.schema)
        writer.write_table(results_table)
        writer.close()
    
    return output_file

# Example usage
if __name__ == "__main__":
    # Define x values for evaluation
    x = np.linspace(20, 40, 100)
    
    # Input and output file paths
    input_file = "/Users/n.t.wamsley/Projects/Pioneer.jl/data/library_spline_tests/example_library_spline_coef.arrow"
    output_file = "/Users/n.t.wamsley/Projects/Pioneer.jl/data/library_spline_tests/evaluated_splines_scipy.arrow"
    
    # Read, evaluate, and save splines
    saved_file = read_evaluate_and_save_splines(input_file, output_file, x)
    print(f"Results saved to: {saved_file}")
    
    # Verify the saved data
    verification_table = pa.ipc.RecordBatchFileReader(
        pa.memory_map(saved_file, 'r')
    ).read_all()
    verification_df = verification_table.to_pandas()
    
    print("\nVerification of saved data:")
    print(f"Number of splines saved: {len(verification_df)}")
    print("\nFirst row:")
    print(f"Annotation: {verification_df.annotation.iloc[0]}")
    print(f"First few spline values: {verification_df.spline_values.iloc[0][:5]}")
    print(f"First few x values: {verification_df.x_values.iloc[0][:5]}")
=#
Arrow.write("/Users/n.t.wamsley/Projects/Pioneer.jl/data/example_library_spline_coef.arrow", df[1:10,:])
df = DataFrame(Arrow.Table("/Users/n.t.wamsley/Projects/Pioneer.jl/data/library_spline_tests/example_library_spline_coef.arrow"))
knots = [6.0, 13.0, 20.0, 27.0, 34.0, 41.0, 48.0, 55.0]
degree = 3
xbins = LinRange(20.0, 40.0, 100)
spline_eval_dict = Dict{String, Vector{Float64}}()
for i in range(1, size(df, 1))
    coefs = Float64.(collect(df[i,:coefficients]))
    spl = BSpline(knots, coefs, 3)
    spline_eval_dict[df[i,:annotation]] = spl.(xbins)
end
scipy_spline_eval = Arrow.Table("/Users/n.t.wamsley/Projects/Pioneer.jl/data/library_spline_tests/evaluated_splines_scipy.arrow")
scipy_eval_dict = Dict{String, Vector{Float64}}()
for i in range(1, length(scipy_spline_eval[:spline_values]))
    scipy_eval_dict[scipy_spline_eval[:annotation][i]] = coalesce.(scipy_spline_eval[:spline_values][i], 0.0)
end

plot(xbins, spline_eval_dict["y3^1"])
plot!(xbins, scipy_eval_dict["y3^1"])
@testset "scipy_splines_compare" begin
    for frag_name in keys(spline_eval_dict)
        @test maximum(
            abs.(
                (spline_eval_dict[frag_name] .- scipy_eval_dict[frag_name])./spline_eval_dict[frag_name]
                )
                ) < 1e-6
    end
end

p = plot()
for frag_name in keys(spline_eval_dict)
    plot!(p, xbins, spline_eval_dict[frag_name], label = nothing, show = true)
    plot!(p, xbins, scipy_eval_dict[frag_name], label = frag_name, show = true)
end

