import sys
import tskit

args = sys.argv
ts = tskit.load(args[1])
output = args[2]

ts_simplified = ts.simplify(keep_unary=False)
ts_simplified.dump(output)