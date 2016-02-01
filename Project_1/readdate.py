from matplotlib import pyplot as plt


myfile = open("energies.txt", "r")
data = []
for line in myfile:
    data.append(float(line));
avarage = sum(data)/len(data)
plt.title(avarage)
plt.plot(data)
plt.axhline(y=avarage, linewidth=2, color='r')
plt.show()
