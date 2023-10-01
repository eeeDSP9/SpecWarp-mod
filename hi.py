import os
import subprocess

# input directory containing audio files to process
input_dir = "C:\\Users\\I5 - 8th GEN\\Desktop\\asr\\input"

# output directory for processed audio files
output_dir = "C:\\Users\\I5 - 8th GEN\\Desktop\\asr\\output"

# path to Matlab executable
matlab_exe = "C:\\Program Files\\MATLAB\\R2022b\\bin\\matlab.exe"

# path to Matlab script for processing audio files
matlab_script = "C:\\Users\\I5 - 8th GEN\\Desktop\\asr\\fm_warp.m"


# iterate over each audio file in the input directory
for filename in os.listdir(input_dir):
    if filename.endswith(".wav"):
        # construct paths to input and output audio files
        input_path = os.path.join(input_dir, filename)
        output_path = os.path.join(output_dir, filename)
        
        # call Matlab script with input audio file and output audio file paths as arguments
        subprocess.run([matlab_exe, "-nodisplay", "-r", f"fm_warp('{input_path}', '{output_path}'); exit;"])