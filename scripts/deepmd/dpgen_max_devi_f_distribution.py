from scipy.stats import gaussian_kde
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 14
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# 已知数据
data = np.loadtxt('./Max_Devi_F')

# 温度
temperatures = [300, 450, 600, 750]

# 绘制所有温度的高斯分布曲线
plt.figure(figsize=(10, 6))

for i, temp in enumerate(temperatures):
    # 使用gaussian_kde估计数据分布
    density = gaussian_kde(data[:, i])
    xs = np.linspace(0, max(data[:, i])*1.1, 200)
    plt.plot(xs, density(xs), label=f'{temp}K')

plt.xlabel("Max_Devi_F Value")
plt.ylabel("Density")
plt.xlim(0,)
plt.ylim(0,)
plt.legend()
plt.grid(True)
plt.savefig('./Max_Devi_F.png', bbox_inches='tight')
