import os 

parentfolder = "/Users/shaeermoeed/Github/DVRPairMC/Results/18_03_2024_12_32_29"

a = os.listdir(parentfolder)
print(a)
print(a[0].split("_"))

a = [(1,2),(3,4),(5,6)]
b = list(map(lambda x :x[1], a))
print(b)
