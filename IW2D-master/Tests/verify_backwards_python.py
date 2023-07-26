# Requires compiled c++ executables in directory ImpedanceWake2D
# Prints outputs to console, since capturing to e.g. a file mixes the order of the outputs from this file and from the actual scripts

import subprocess
import os
from pathlib import Path
from filecmp import cmpfiles

# import sys
# sys.stdout = open("test_results.txt", "w")

examples_path = Path(__file__).resolve().parents[1] / "examples"

tests = [
    (examples_path / "input_files/RoundChamberInputFile.txt", "roundchamber"),  # test roundchamber
    (examples_path / "input_files/FlatChamberInputFile.txt", "flatchamber"),    # test flatchamber with DC resistivity + relaxation time and susceptebility + relaxation time
    (examples_path / "input_files/test_empty.txt", "flatchamber"),              # test a completely empty file
    (examples_path / "input_files/test_tandelta.txt", "flatchamber"),           # test tandelta and multiple layers
    (examples_path / "input_files/test_freq_eps_and_mu.txt", "flatchamber"),    # test freq dependent complex eps1 and mu1
    (examples_path / "input_files/test_freq_sigma.txt", "flatchamber"),         # test freq dependent complex sigma
    (examples_path / "input_files/FlatChamberInputFileNonlinear.txt", "flatchamber"),   # test higher order flatchamber
    (examples_path / "input_files/test_lin_round.txt", "roundchamber"),         # test linear frequency scan, and provoke strange result mismatch
]

total_matches = 0
total_mismatches = 0
total_errors = 0
total_line_mismatches = 0

# set working directory and create folder for outputs
os.chdir(Path(__file__).resolve().parents[1])
if "output_cpp" not in os.listdir():
    os.mkdir("output_cpp")
if "output_py" not in os.listdir():
    os.mkdir("output_py")



for inputfile, mode in tests:

    print(f"Running c++ {mode} executable with input file {inputfile}")
    with open(inputfile) as f:
        subprocess.run(f"./IW2D/cpp/{mode}.x", stdin=f)
    for file in os.listdir("."):
        if file.endswith(".dat"):
            os.replace(file, "./output_cpp/"+file)

    print()

    print(f"Running python module in {mode} mode with input file {inputfile}")
    subprocess.run(f"python3 -m IW2D.legacy {mode} {inputfile}".split())
    for file in os.listdir("."):
        if file.endswith(".dat"):
            os.replace(file, "./output_py/"+file)
    
    print()

    match, mismatch, error = cmpfiles("output_cpp", "output_py", os.listdir("output_py"))
    print(f"{len(match)} matching files:")
    for m in match:
        print("\t",m)
        total_matches += 1
    print(f"{len(mismatch)} non-matching files:")
    for mm in mismatch:
        print("\t",mm)
        total_mismatches += 1
    print(f"{len(error)} errors:")
    for e in error:
        print("\t",e)
        total_errors  += 1
    
    print()

    for mm in mismatch:
        with open("output_cpp/"+mm, "r") as fcpp:
            with open("output_py/"+mm, "r") as fpy:
                for i, (lcpp, lpy) in enumerate(zip(fcpp.readlines(), fpy.readlines())):
                    if lcpp != lpy:
                        total_line_mismatches += 1
                        print(f"Difference in line {i} of {mm}")
                        print("cpp result: ", lcpp.strip())
                        print("py result:  ", lpy.strip())
                        print()
    
    print(("="*116) + "\n")

    for fil in os.listdir("output_cpp"):
        os.remove("output_cpp/"+fil)
    for fil in os.listdir("output_py"):
        os.remove("output_py/"+fil)

print(f"""\
Summary:
    matching files:     {total_matches}
    mismatching files:  {total_mismatches}
    errors:             {total_errors}
    total lines diff:   {total_line_mismatches}\t(note: the same frequency will often contribute more than once here)\
""")

os.rmdir("output_cpp")
os.rmdir("output_py")

# sys.stdout.close()
