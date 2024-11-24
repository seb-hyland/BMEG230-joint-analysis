module JointAnalysis

using CSV
using DataFrames
using LinearAlgebra
using Plots
using SavitzkyGolay


#=
Calculate joint angle over time from 3 points, in 3 dimensions

### Parameters:
- `a_x`, `a_y`, `a_z`: Vectors representing the coordinates of the first point
- `b_x`, `b_y`, `b_z`: Vectors representing the coordinates of the second point
- `c_x`, `c_y`, `c_z`: Vectors representing the coordinates of the third point

### Returns:
- A `Vector{Float64}` describing the joint angle over time

### Requirements:
- All input vectors have the same length

=#
function calculate_angles(a_x, a_y, a_z, b_x, b_y, b_z, c_x, c_y, c_z)
    angles = Float64[]

    for i in 1:length(a_x)
        vec_1 = [b_x[i] - a_x[i], b_y[i] - a_y[i], b_z[i] - a_z[i]]
        vec_2 = [c_x[i] - b_x[i], c_y[i] - b_y[i], c_z[i] - b_z[i]]
        
        norm_1 = normalize(vec_1)
        norm_2 = normalize(vec_2)
        
        dot_product = dot(norm_1, norm_2)
        angle = acos(dot_product)
        push!(angles, rad2deg(angle))
    end

    return angles
end


#=
Process data to plot knee, hip, and elbow angles over time, and plot the data with respect to time

### Parameters:
- `dir::String`:        The path to the data directory
- `file_name::String`:  The name of the data file

### Requires:
- `dir` and `file_name` must describe a CSV file
- The CSV file must include a column `frame`
- The CSV file must include columns for all combinations of:
  - `left_JOINT_DIMENSION`, where:
    - `JOINT` is one of `shoulder`, `elbow`, `wrist`, `hip`, `knee`, or `ankle`.
    - `DIMENSION` is one of `x`, `y`, or `z`.

### Outputs:
- Saves a plot in `.svg` format to `./figs/`, named based on the `file_name`.
- Plots assume 30 frames per second

=#
function process_data(dir, file_name)
    # Read CSV data into DataFrame
    file_path = joinpath(dir, file_name)
    data = CSV.read(file_path, DataFrame)

    # Extract relevant data
    frame = data.frame

    S_x = data.left_shoulder_x
    S_y = data.left_shoulder_y
    S_z = data.left_shoulder_z

    E_x = data.left_elbow_x
    E_y = data.left_elbow_y
    E_z = data.left_elbow_z

    W_x = data.left_wrist_x
    W_y = data.left_wrist_y
    W_z = data.left_wrist_z

    H_x = data.left_hip_x
    H_y = data.left_hip_y
    H_z = data.left_hip_z

    K_x = data.left_knee_x
    K_y = data.left_knee_y
    K_z = data.left_knee_z

    A_x = data.left_ankle_x
    A_y = data.left_ankle_y
    A_z = data.left_ankle_z

    # Calculate angles over time
    knee_angles = calculate_angles(H_x, H_y, H_z, K_x, K_y, K_z, A_x, A_y, A_z)
    hip_angles = calculate_angles(S_x, S_y, S_z, H_x, H_y, H_z, K_x, K_y, K_z)
    elbow_angles = calculate_angles(S_x, S_y, S_z, E_x, E_y, E_z, W_x, W_y, W_z)

    # Smooth data
    sm_knee_angles = savitzky_golay(knee_angles, 9, 3).y
    sm_hip_angles = savitzky_golay(hip_angles, 9, 3).y
    sm_elbow_angles = savitzky_golay(elbow_angles, 9, 3).y

    # Convert frames to time
    time = frame / 30
    
    # Plot
    p = plot(time, [sm_knee_angles, sm_hip_angles, sm_elbow_angles],
             xlabel = "Time (Frame)",
             ylabel = "Joint Angle (Degrees)",
             title = "Joint Angles Over Time",
             fontfamily = "Courier",
             linewidth = 2,
             thickness_scaling = 1,
             palette = ["#FF4747", "#00BCB4", "#FFB547"],
             label = ["Knee Angle" "Hip Angle" "Elbow Angle"])
    savefig(p, "./figs/$(splitext(file_name)[1])_plot.svg")
end


#= -----------------------------------------------------------------------------------
Main script, which iterates over every CSV file in the data directory and processes it
----------------------------------------------------------------------------------- =#
directory = "./data/"
files = readdir(directory)
csv_files = filter(f -> endswith(f, ".csv"), files)
if !isdir("./figs")
    mkdir("./figs")
end
for csv_file in csv_files
    process_data(directory, csv_file)
end

end
