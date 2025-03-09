#!/bin/bash

# Input parameters file
input_file="parameters_MMS.txt"

sed -i "s/^C_e = .*/C_e = 20                                        /" "$input_file"
sed -i "s/^D_eff = .*/D_eff = 1e-10                                 /" "$input_file"
sed -i "s/^R = .*/R = 0.5                                           /" "$input_file"
sed -i "s/^N = .*/N = 1000                                          /" "$input_file"
sed -i "s/^k = .*/k = 4e-9                                          /" "$input_file"
sed -i "s/^T = .*/T = 4e9                                           /" "$input_file"
sed -i "s/^K = .*/K = 10000                                         /" "$input_file"
sed -i "s/^C_i = .*/C_i = 20                                        /" "$input_file"
sed -i "s/^source = .*/source = 1                                   /" "$input_file"
sed -i "s/^output_filename = .*/output_filename = output_MMS.csv    /" "$input_file"

########################################## Space convergence ##########################################
# Modify the K value in the parameters file
sed -i "s/^K = .*/K = 10000        # Nombre de pas de temps [-]/" "$input_file"

# Output file to store convergence results
output_file="convergence_results_space.csv"

# Range of N values to test
N_values=(10 20 50 100 200 500 1000 2000 5000 10000)

# Header for the output file
echo "N,L1_error,L2_error,L_inf_error" > "$output_file"

# Loop over N values
for N in "${N_values[@]}"; do
    echo "Running simulation with N = $N"

    # Modify the N value in the parameters file
    sed -i "s/^N = .*/N = $N        # Nombre de noeuds [-]/" "$input_file"

    # Run the C++ program and capture the output
    output=$(./main "$input_file")

    # Extract errors from the output
    L1_error=$(echo "$output" | grep "L1 error" | awk '{print $4}')
    L2_error=$(echo "$output" | grep "L2 error" | awk '{print $4}')
    L_inf_error=$(echo "$output" | grep "L_inf error" | awk '{print $4}')

    # Append results to the output file
    echo "$N,$L1_error,$L2_error,$L_inf_error" >> "$output_file"
done

echo "Convergence analysis complete. Results saved to $output_file"

########################################## Time convergence ##########################################
# Modify the N value in the parameters file
sed -i "s/^N = .*/N = 1000        # Nombre de noeuds [-]/" "$input_file"

# Output file to store convergence results
output_file="convergence_results_time.csv"

# Range of N values to test
K_values=(10 20 50 100 200 500 1000 2000 5000 10000)

# Header for the output file
echo "K,L1_error,L2_error,L_inf_error" > "$output_file"

# Loop over N values
for K in "${K_values[@]}"; do
    echo "Running simulation with K = $K"

    # Modify the K value in the parameters file
    sed -i "s/^K = .*/K = $K        # Nombre de pas de temps [-]/" "$input_file"

    # Run the C++ program and capture the output
    output=$(./main "$input_file")

    # Extract errors from the output
    L1_error=$(echo "$output" | grep "L1 error" | awk '{print $4}')
    L2_error=$(echo "$output" | grep "L2 error" | awk '{print $4}')
    L_inf_error=$(echo "$output" | grep "L_inf error" | awk '{print $4}')

    # Append results to the output file
    echo "$K,$L1_error,$L2_error,$L_inf_error" >> "$output_file"
done