import sys
import tsinfer

args = sys.argv
ancestors = tsinfer.load(args[1])
#uprtime = args[2]
#lwertime = args[3]
output = args[2]

ancestors.truncate_ancestors(
    path=output, 
    max_file_size=2**38,
    lower_time_bound=0.2,
    upper_time_bound=0.4,
    length_multiplier=2,
)


