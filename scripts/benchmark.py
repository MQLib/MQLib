# The following script performs the LINPACK 100 benchmark to compute a flop
# count.
# 
# To run the LINPACK 100 benchmark on an EC2 compute node, do the following:
# 1) Start up your instance and log in
# 2) Copy benchmark.py to the compute node
# 3) Install fortran:
#
#      # Linux:
#      > sudo apt-get install gfortran
#      # Mac:
#      > brew install gcc
# 
# 4) Download the LINPACK 100 fortran code:
# 
#      > wget http://www.netlib.org/benchmark/linpackd -O linpackd.f
# 
# 5) Compile the LINPACK 100 fortran code with maximum compiler optimization:
# 
#      > gfortran linpackd.f -O3
# 
# 6) Run the benchmark repeatedly and get the flop count by averaging the runtime
# 
#      > python benchmark.py
# 
# 7) Now you have benchmark results. The first and second numbers are based on
#    "times for array with leading dimension 201" and "times for array with
#    leading dimension 200", respectively. According to
#    http://www.netlib.org/utk/people/JackDongarra/faq-linpack.html, question
#    "How can I interpret the results from the Linpack 100x100 benchmark?", 
#    "This is done to see what effect, if any, the placement of the arrays in
#    memory has on the performance." I think we should use the faster one, which
#    appears to be the second. Here are the results on m3.medium, run with
#    10,000 benchmarks:
# 
#       Flop count 1 (Mflop): 955.466984757
#       Flop count 2 (Mflop): 1128.95224688
# 

import subprocess

all1 = []
all2 = []
ops = 2.0/3.0*100.0**3 + 2.0*100.0**2
for rep in range(10000):
    proc = subprocess.Popen(["./a.out"], stdout=subprocess.PIPE)
    output = proc.stdout.read().decode("utf-8")
    all1 += [float(output.split("\n")[x].strip().split()[2]) for x in [20, 21, 22, 23]]
    all2 += [float(output.split("\n")[x].strip().split()[2]) for x in [26, 27, 28, 29]]

print("Flop count 1 (Mflop):", ops * len(all1) / sum(all1) / 1e6)
print("Flop count 2 (Mflop):", ops * len(all2) / sum(all2) / 1e6)
