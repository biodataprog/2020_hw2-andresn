#!/usr/bin/env python3
#%%
start = 0
end = 99
divisor = 7
print("Printing out numbers from", start, "to", end, " not divisible by", divisor)
print('\n'.join(str(x) for x in range(start, end+1) if x % divisor))
#%%