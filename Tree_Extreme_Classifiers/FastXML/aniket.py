import numpy as np

ret = np.zeros( 3400 )
i = -1
f = open("train.y", "r")
for temp in f:
    if i==-1:
        i += 1
        continue
    lol = [int(y[0]) for y in  [x.split(":") for x in str(temp).split(" ")]]
    for x in lol:
        ret[x] += 1
f.close()
total = sum(ret) 
ret = [(1 + 14.368*((x + 1.5)**(-0.55))) for x in ret]
for t in ret:
	print(t);