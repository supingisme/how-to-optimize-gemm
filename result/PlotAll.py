#!/usr/bin/env python
import matplotlib.pyplot as plt

import numpy as np

# Indicate the number of floating point operations that can be executed
# per clock cycle
nflops_per_cycle = 4

# Indicate the number of processors being used (in case you are using a
# multicore or SMP)
nprocessors = 1

# Indicate the clock speed of the processor.  On a Linux machine this info
# can be found in the file /proc/cpuinfo
#
# Note: some processors have a "turbo boost" mode, which increases
# the peak clock rate...
#
GHz_of_processor = 2.419


class Parser:
    def __init__(self, file_name=None):
        self.attrs = {}
        with open(file_name) as file:
            self.toks = file.read().split()
            self.toksi = 0
            file.close()
            self.attrs = self.parse()

    def next(self):
        tok = self.toks[self.toksi]
        self.toksi += 1
        return tok

    def get_var_name(self):
        return self.next()

    def get_symbol(self, sym):
        tok = self.next()
        assert(tok == sym)
        return tok

    def get_value(self):
        value = None
        tok = self.next()
        if tok == '[':
            # list
            value = []
            tok = self.next()
            while not tok.startswith(']'):
                value.append(float(tok))
                tok = self.next()
        elif tok.startswith("'"):
            value = tok[1:-2]
        
        assert value != None
        return value
    
    def parse(self):
        res = {}
        while self.toksi < len(self.toks):
            var = self.get_var_name()
            self.get_symbol('=')
            val = self.get_value()
            res[var] = val
        return res
    
    def __getattr__(self, name):
        return self.attrs[name]

test1 = Parser("test_MMult0.m")
test2 = Parser("test_MMult1.m")
test3 = Parser("test_MMult_1x4_3.m")
test4 = Parser("test_MMult_1x4_4.m")
test5 = Parser("test_MMult_1x4_5.m")
test6 = Parser("test_MMult_1x4_6.m")
test7 = Parser("test_MMult_1x4_7.m")
test8 = Parser("test_MMult_1x4_8.m")
test9 = Parser("test_MMult_1x4_9.m")
test10 = Parser("test_MMult_4x4_9.m")
test11 = Parser("test_MMult_4x4_11.m")
test12 = Parser("test_MMult_4x4_12.m")
test13 = Parser("test_MMult_4x4_13.m")
# test14 = Parser("test_MMult_4x4_14.m")
test15 = Parser("test_MMult_4x4_15.m")

# print(test1)

# print(test1.MY_MMult);
test1_data = np.array(test1.MY_MMult).reshape(-1, 6)
test2_data = np.array(test2.MY_MMult).reshape(-1, 6)
test3_data = np.array(test3.MY_MMult).reshape(-1, 6)
test4_data = np.array(test4.MY_MMult).reshape(-1, 6)
test5_data = np.array(test5.MY_MMult).reshape(-1, 6)
test6_data = np.array(test6.MY_MMult).reshape(-1, 6)
test7_data = np.array(test7.MY_MMult).reshape(-1, 6)
test8_data = np.array(test8.MY_MMult).reshape(-1, 6)
test9_data = np.array(test9.MY_MMult).reshape(-1, 6)
test10_data = np.array(test10.MY_MMult).reshape(-1, 6)
test11_data = np.array(test11.MY_MMult).reshape(-1, 6)
test12_data = np.array(test12.MY_MMult).reshape(-1, 6)
test13_data = np.array(test13.MY_MMult).reshape(-1, 6)
# test14_data = np.array(test14.MY_MMult).reshape(-1, 6)
test15_data = np.array(test15.MY_MMult).reshape(-1, 6)
# print(test1_data)

max_gflops = nflops_per_cycle * nprocessors * GHz_of_processor;

fig, ax = plt.subplots()
ax.plot(test1_data[:,0], test1_data[:,1], 'r', label=test1.version)
ax.plot(test2_data[:,0], test2_data[:,1], 'g', label=test2.version)
ax.plot(test3_data[:,0], test3_data[:,1], 'b', label=test3.version)
ax.plot(test4_data[:,0], test4_data[:,1], 'c', label=test4.version)
ax.plot(test5_data[:,0], test5_data[:,1], 'm', label=test5.version)
ax.plot(test6_data[:,0], test6_data[:,1], 'y', label=test6.version)
ax.plot(test7_data[:,0], test7_data[:,1], 'k', label=test7.version)
ax.plot(test8_data[:,0], test8_data[:,1], 'ro-', label=test8.version)
ax.plot(test9_data[:,0], test9_data[:,1], 'go-', label=test9.version)
ax.plot(test10_data[:,0], test10_data[:,1], 'bo-', label=test10.version)
ax.plot(test11_data[:,0], test11_data[:,1], 'co-', label=test11.version)
ax.plot(test12_data[:,0], test12_data[:,1], 'mo-', label=test12.version)
ax.plot(test13_data[:,0], test13_data[:,1], 'yo-', label=test13.version)
ax.plot(test15_data[:,0], test15_data[:,1], 'ko-', label=test15.version)

ax.set(xlabel='m = n = k', ylabel='GFLOPS/sec.',
       title="Calculate ability")
ax.grid()
ax.legend()

ax.set_xlim([test1_data[0,0], test2_data[-1,0]])
ax.set_ylim([0, max_gflops])

manager = plt.get_current_fig_manager()
manager.resize(*manager.window.maxsize())

plt.show()
fig.savefig("result.png")
