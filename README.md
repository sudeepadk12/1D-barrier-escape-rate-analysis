# 1D-barrier-escape-rate-analysis
# Description
This repository contains one dimensional analysis of thermally activated biased barrier crossing analysis in biopolymers. The analysis gives critical force/ rupture force for unfolding of biopolymers when time-depending pulling force is applied. The unfolding of biopolymer is equivalent to the extinction of the barrier.

# Features
- Simulates multiple pulling rates.
  - Supports custom pulling rates
  - Handles multiple rates in single simulation
- Runs multiple (```5``` in this example) iterations for each pulling rate.
  - Independent random initialization for each run
  - Consistent output formatting for compatibility
  
- Outputs results in a file (```results.dat```) with the format:
  - Run Number
  - Pulling Rate
  - External Force
# Prerequisites
Ensure you have installed the following
- GCC or any C++17 compatible compiler
  - Recommended: GCC 9.0 or higher for optimal performance.
- Make (optional but recommended for building)
  - Ensures easier compilation and artifact management.
 
# Files
- ```modver.cpp```: Main Source Code
- ```Makefile```: Instructions for building and running the simulation
  - Automates the compilation and clean up
  - Includes a ```run``` target for direct execution
- ```results.dat```: Output file for simulation result
  - Tabular format for easier critical force analysis and plotting
# Compilation and Execution
To compile and run the program
1. Build the program
  - Command:   
  ```make```
  - This compiles the source code (```modver.cpp```) and generates an executable
3. Run the simulation
  - Command:    
  ``` make run```
  - This executes the program and generates the output file (```results.dat```)
4. Clean build artifacts
  - Command:   
  ```make clean ```
  - This removes the executable and other temporary files, ensuring a clean working directory.
# Output
The results are saved in the ```results.dat``` file. Each row contains:
1. The adjustable run numbers ( 1 to 5 in this example)
   - Minimum 2500 runs required for a smooth probabilioty distribution or cumulative probability distribution of critical force for a particular pulling rate
2. The value of pulling rates ( e.g., 0.0001, 0.001, etc.).
   - Corresponds to the pulling rate at which external force is applied.
3. The computed force.
   - Critical force at which the particle escapes the barrier or the biopolymer unfolds.
   
# Example Output 
```Run_Num```   ```Pulling_Rate```  ```External_Force```  
```1 ```              ``` 0.0001 ```              ```  12.67 ```   
```2 ```              ``` 0.0001 ```              ```  14.35 ```  
``` ...```


# Notes
- Modify the pulling_rates vector in the source code to add or remove pulling rates as needed.
  - Example: Add additional rates like 0.00001, 0.00002 for slower pulling.
- The PRECISION and other constants can be adjusted in the source code for finer control of the simulation.
  - Smaller precision values may improve accuracy but increase runtime.
# License
This project is open-source. Use, modify, and distribute freely.

