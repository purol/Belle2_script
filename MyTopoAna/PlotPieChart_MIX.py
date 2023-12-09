import matplotlib.pyplot as plt
import numpy as np

# [('others', 8), ('511=>-413|211|211|-211|111|', 6), ('511=>213|-411|', 3), ('511=>213|-413|', 2), ('511=>-413|20213|', 2), ('511=>-411|221|211|', 1), ('511=>-411|211|-213|211|22|', 1), ('511=>-411|223|211|', 1), ('511=>-411|113|211|111|', 1), ('511=>213|-411|22|', 1), ('511=>-411|-211|211|211|', 1), ('511=>-413|321|-22|', 1), ('511=>20213|-411|', 1)]

def autopct_format(values):
    def my_format(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{v:d}'.format(v=val)
    return my_format

Category = [
r"$B^{0}\rightarrow D^{*}(2010)^{-} \mu^{+} \nu_{\mu}$", # 511=>-413|-13|14| 130
r"$B^{0}\rightarrow D^{-} \mu^{+} \nu_{\mu}$", # 511=>-411|-13|14| 130
r"$B^{0}\rightarrow D^{*}(2010)^{-} e^{+} \nu_{e}$", # 511=>-413|-11|12| 109
r"$B^{0}\rightarrow D^{-} e^{+} \nu_{e}$", # 511=>-411|-11|12| 108
r"$B^{0}\rightarrow D^{*}(2010)^{-} a_{1}(1260)^{+}$", # 511=>-413|20213| 50
r"$B^{0}\rightarrow D^{*}(2010)^{-} \pi^{+} \pi^{0}$", # 511=>-413|211|111| 48
r"$B^{0}\rightarrow \rho(770)^{+} D^{-}$", # 511=>213|-411| 47
r"$B^{0}\rightarrow \rho(770)^{+} D^{*}(2010)^{-}$", # 511=>213|-413| 47
r"$B^{0}\rightarrow D^{*}(2010)^{-} \pi^{+}$", # 511=>-413|211| 45
r"$B^{0}\rightarrow D^{*}(2010)^{-} \pi^{+} \pi^{+} \pi^{-} \pi^{0}$", # 511=>-413|211|211|-211|111| 37
r"$B^{0}\rightarrow D^{-} \tau^{+} \nu_{\tau}$", # 511=>-411|-15|16| 35
r"$B^{0}\rightarrow D^{*}(2010)^{-} \rho(770)^{0} \pi^{+}$", # 511=>-413|113|211| 35
r"$B^{0}\rightarrow D^{*}(2010)^{-} \tau^{+} \nu_{\tau}$", # 511=>-413|-15|16| 33
"others"]

TotalN = 1622
values = [130, 130, 109, 108, 50, 48, 47, 47, 45, 37, 35, 35, 33]
values.append(TotalN - sum(values))

plt.rcParams['figure.figsize'] = [10,10]
plt.rcParams["font.size"] = 12.0

cmap = plt.get_cmap('Spectral')
colors = [cmap(i) for i in np.linspace(0,1,len(values))]
plt.pie(values, labels=Category, colors = colors, autopct=autopct_format(values), radius = 0.8)

plt.savefig("MIX_BKG_PIE.png")
