# Writes js/results.js for the web page from the two JSON result files
# (run ex1_integration.py and ex2_normal_cdf.py first).
import json, mpmath
mpmath.mp.dps = 40
e1 = json.load(open("ex1_results.json")); e2 = json.load(open("ex2_results.json"))
ref = [mpmath.nstr(mpmath.ncdf(mpmath.mpf(i) / 100), 20) for i in range(-1000, 1001)]
with open("../js/results.js", "w") as f:
    f.write("// numbers produced by code/ex1_integration.py and code/ex2_normal_cdf.py\n")
    f.write("const EX1 = " + json.dumps(e1) + ";\n")
    f.write("const EX2 = " + json.dumps(e2) + ";\n")
    f.write("// reference values of Phi(z) for z = -10.00, -9.99, ..., 10.00 (mpmath, 40 digits, rounded to 20)\n")
    f.write("const PHI_REF = [" + ",".join(ref) + "];\n")
