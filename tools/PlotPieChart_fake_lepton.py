import matplotlib.pyplot as plt
import numpy as np

# [('521=>213|-421|', 11), ('others', 7), ('521=>-423|211|', 6), ('521=>-413|111|211|211|', 6), ('521=>-423|20213|', 6), ('521=>-423|-211|211|211|111|', 4), ('521=>-421|211|', 4), ('521=>-423|213|', 3), ('521=>-413|213|211|', 2), ('521=>-421|113|211|111|', 2), ('521=>-421|113|211|', 2), ('521=>-413|211|211|', 2), ('521=>-411|111|211|211|', 2), ('521=>-421|-211|211|211|', 2), ('521=>-423|211|-211|211|', 2), ('521=>-423|113|211|111|', 1), ('521=>-421|223|211|', 1), ('521=>-411|213|211|', 1), ('521=>20213|-421|', 1), ('521=>-421|113|211|22|', 1), ('521=>-421|431|', 1), ('521=>213|211|-211|211|-211|', 1), ('521=>321|-213|-443|211|', 1), ('521=>-423|-223|211|', 1), ('521=>-423|-221|211|113|', 1), ('521=>-423|321|', 1), ('521=>-423|-211|211|211|111|22|', 1)]

def autopct_format(values):
    def my_format(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{v:d}'.format(v=val)
    return my_format

Category = [
r"$K$ from $K$",
r"$K$ from $\pi$",
r"$K$ from $e$", 
r"$K$ from $\mu$"]

TotalN = 238.166698
values = [238.166698, 29.878096, 0.163917, 0.0]
#values.append(TotalN - sum(values))

plt.rcParams['figure.figsize'] = [10,10]
plt.rcParams["font.size"] = 12.0

#cmap = plt.get_cmap('Spectral')
#colors = ['lightblue', 'deepskyblue', 'turquoise', 'cornsilk', 'tan']
#plt.pie(values, labels=Category, explode = [0.0, 0.2, 0.35, 0.48, 0.59], colors = colors, autopct='%.1f%%', radius = 0.8, startangle = 30)

cmap = plt.get_cmap('Spectral')
colors = ['lightblue', 'deepskyblue', 'turquoise', 'cornsilk']
plt.pie(values, labels=Category, explode = [0.0, 0.2, 0.35, 0.48], colors = colors, autopct='%.1f%%', radius = 0.8, startangle = 30)


plt.savefig("fake_PIE.png")
