#!/bin/bash

# Input parameters file
input_file="parameters.txt"

# Modify the K value in the parameters file
sed -i "s/^K = .*/K = 10        # Nombre de pas de temps [-]/" "$input_file"

# Output file to store convergence results
output_file="convergence_results.csv"

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