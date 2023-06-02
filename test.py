import time
import timeit
import statistics as stat
import matplotlib.pyplot as plt
import numpy as np
import sys
from io import StringIO

#usage: python test.py <number of repeats>
#           -i <module_1> .. <module_n>
#           -f <function> {OR} -f <function_1> .. <function_n>
#           -a <arg_1> .. <arg_m>
#           -e <extra modules here> (optional)

#parse command line args
reps = int(sys.argv[1])
modules = []
funcs = []
args = []
extra = []

j = 2
argv_len = len(sys.argv)
while j < argv_len:
    if sys.argv[j] == '-i':
        j += 1
        while j < argv_len and not sys.argv[j][0] == '-':
            modules.append(sys.argv[j])
            j += 1
    elif sys.argv[j] == '-f':
        j += 1
        while j < argv_len and not sys.argv[j][0] == '-':
            funcs.append(sys.argv[j])
            j += 1
    elif sys.argv[j] == '-a':
        j += 1
        while j < argv_len and not sys.argv[j][0] == '-':
            args.append(sys.argv[j])
            j += 1
    elif sys.argv[j] == '-e':
        j += 1
        while j < argv_len and not sys.argv[j][0] == '-':
            extra.append(sys.argv[j])
            j += 1
n = len(modules)
if (len(funcs) == 1):
    for i in range(1, n):
        funcs.append(funcs[0])

#make tests
test = []
setup = []
data = []

for i in range(n):
    s = '''\
from {module} import {func} as func
    '''.format(module=modules[i], func=funcs[i])
    for i in range(len(extra)):
        s += '\nimport ' + extra[i] + '\n'
    

    t = 'func('
    for i in range(len(args)):
        t += args[i]
        if not i == len(args)-1:
            t += ','
        else:
            t += ')'
    #print(s)
    #print(t)

    setup.append(s)
    test.append(t)

for i in range(n):
    l = timeit.repeat(test[i], setup[i],
                      repeat = reps,
                      number = 1)
    data.append(l)

#calculate minimum/median of data
bar1 = []
bar2 = []
for i in range(n):
    bar1.append(min(data[i]))

for i in range(n):
    bar2.append(stat.median(data[i]))

#configure horizontal bar plot
fig, ax = plt.subplots()
ind = np.arange(n)
width = 0.35
p1=ax.barh(ind,bar1,label='min',height=0.35)
p2=ax.barh(ind+width,bar2,label='median', height=0.35)

ax.set_xlabel('seconds')
ax.set_yticks(ind)
ax.set_yticklabels(modules)
ax.legend(loc = 'lower right')
ax.invert_yaxis()

#make title
title = 'Execution time of '
for i in range(n):
    title += modules[i]
    if not i == n - 1:
        title += ' vs. '
    else:
        title += ' in seconds'
ax.set_title(title)

#capture output for comparison
compare = []
for i in range(n):
    tmp = sys.stdout
    my_result = StringIO()
    sys.stdout = my_result
    code = setup[i] + '\nprint(' + test[i] + ')'
    exec(code)
    sys.stdout = tmp
    compare.append(my_result.getvalue())

#print output
print('Output:')
print('1: '+ modules[0] + ' (reference implementation):')
print('\t' + compare[0], end='')

for i in range(1, n):
    print(str(i+1) + ': ' + modules[i] + ':')
    print('\t' + compare[i], end='')

#show plot
plt.show()
