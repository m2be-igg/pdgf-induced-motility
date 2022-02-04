#!/bin/bash

# Define number of replicates
numberOfReplicates=3

# Make a directory to run simulations and another to store the results
mkdir -p final_output

# Loop through the replicates
for i in $(seq $numberOfReplicates)
do
   # (change project3D to the name of the PhysiCell project you want to run)
   project.exe

   # Store the results in the general results folder
   mkdir ./final_output/output$i
   cp -r output/. ./final_output/output$i

   # Clean output folder to avoid mixing data if something fails
   rm -r output
   mkdir output

done