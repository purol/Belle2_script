import matplotlib.pyplot as plt
import numpy as np

# [('others', 8), ('511=>-413|211|211|-211|111|', 6), ('511=>213|-411|', 3), ('511=>213|-413|', 2), ('511=>-413|20213|', 2), ('511=>-411|221|211|', 1), ('511=>-411|211|-213|211|22|', 1), ('511=>-411|223|211|', 1), ('511=>-411|113|211|111|', 1), ('511=>213|-411|22|', 1), ('511=>-411|-211|211|211|', 1), ('511=>-413|321|-22|', 1), ('511=>20213|-411|', 1)]

def autopct_format(values):
    def my_format(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{v:d}'.format(v=val)
    return my_format

Category = [r"$B^{0}\rightarrow D^{*}(2010)^{-} \pi^{+} \pi^{+} \pi^{-} \pi^{0}$", # 511=>-413|211|211|-211|111|', 6
r"$B^{0}\rightarrow \rho(770)^{+} D^{-}$", # 511=>213|-411|', 3
r"$B^{0}\rightarrow \rho(770)^{+} D^{*}(2010)^{-}$", # 511=>213|-413|', 2
r"$B^{0}\rightarrow D^{*}(2010)^{-} a_{1}(1260)^{+}$", # 511=>-413|20213|', 2
r"$B^{0}\rightarrow D^{-} \eta \pi^{+}$", # 511=>-411|221|211|', 1
r"$B^{0}\rightarrow D^{-} \pi^{+} \rho(770)^{-} \pi^{+} \gamma$", # 511=>-411|211|-213|211|22|', 1
r"$B^{0}\rightarrow D^{-} \omega(782) \pi^{+}$", # 511=>-411|223|211|', 1
r"$B^{0}\rightarrow D^{-} \rho(770)^{0} \pi^{+} \pi^{0}$", # 511=>-411|113|211|111|', 1
r"$B^{0}\rightarrow \rho(770)^{+} D^{-} \gamma$", # 511=>213|-411|22|', 1
"others"]

TotalN = 29
values = [6, 3, 2, 2, 1, 1, 1, 1, 1]
values.append(TotalN - sum(values))

plt.rcParams['figure.figsize'] = [10,10]
plt.rcParams["font.size"] = 12.0

cmap = plt.get_cmap('Spectral')
colors = [cmap(i) for i in np.linspace(0,1,len(values))]
plt.pie(values, labels=Category, colors = colors, autopct=autopct_format(values), radius = 0.8)

plt.savefig("MIX_BKG_PIE.png")
