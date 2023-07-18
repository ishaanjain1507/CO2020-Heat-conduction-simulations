import csv

# read flash.dat to a list of lists
datContent = [i.strip().split() for i in open("./T_x_0000.dat").readlines()]

# write it as a new CSV file
with open("./T_x_0000.csv", "wb") as f:
    writer = csv.writer(f)
    writer.writerows(datContent)
