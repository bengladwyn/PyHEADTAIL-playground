"""Call signature on command line: "python3 -m IW2D.legacy <mode> <input-file-path>", 
where mode must be one of "flatchamber", "roundchamber", "wake_flatchamber", or "wake_roundchamber"."""

import sys
import cppyy
from time import time

from .flatchamber import flatchamber
from .roundchamber import roundchamber
from .wake_flatchamber import wake_flatchamber
from .wake_roundchamber import wake_roundchamber

VALID_MODES = ["flatchamber", "roundchamber", "wake_flatchamber", "wake_roundchamber"]

# Get command line arguments
mode, input_file = sys.argv[1:]

# Check that supplied mode argument is valid
if mode not in VALID_MODES:
    msg = f"""First argument should be calculation mode.
    Got: '{mode}'
    Expected one of: {VALID_MODES}"""
    raise ValueError(msg)

print(f"Running IW2D in {mode} mode with legacy input format.")
starttime = time()

# Try to open file and put inputs into dictionary
with open(input_file, "r") as f:
    input_dict = {}
    for i, line in enumerate(f.readlines()):
        if not line.isspace():
            if "\t" not in line:
                msg = f"Error reading {input_file}:\n\tNo \\t found on line {i+1}."
                raise ValueError(msg)

            entry_name, value = line.split("\t")
            input_dict[entry_name.strip(":")] = value.strip()

# Do final mode specific inclusions and run main script
if mode == "flatchamber":
    flatchamber(input_dict, input_file)

elif mode == "roundchamber":
    roundchamber(input_dict, input_file)

elif mode == "wake_flatchamber":
    wake_flatchamber(input_dict, input_file)

elif mode == "wake_roundchamber":
    wake_roundchamber(input_dict, input_file)

endtime = time()
print(f"Elapsed time during calculation: {endtime - starttime:.2f} seconds")
