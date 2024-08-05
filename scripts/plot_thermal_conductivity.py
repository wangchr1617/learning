# 文件名: plot_thermal_conductivity.py
# 运行方法: python plot_thermal_conductivity.py
# 功能描述: 加载热导率数据并绘制热导率随时间变化的图。

import numpy as np
import matplotlib.pyplot as plt
from gpyumd.load import load_kappa
from gpyumd.math import running_ave

# 设置matplotlib全局参数
aw = 2
fs = 16
font = {'size': fs}
plt.rc('font', **font)
plt.rc('axes', linewidth=aw)

def set_fig_properties(ax_list):
    """
    设置图形属性。

    参数:
        ax_list (list): Axes对象列表。
    """
    tl = 8  # 大刻度长度
    tw = 2  # 大刻度宽度
    tlm = 4  # 小刻度长度

    for ax in ax_list:
        ax.tick_params(which='major', length=tl, width=tw)
        ax.tick_params(which='minor', length=tlm, width=tw)
        ax.tick_params(which='both', axis='both', direction='in', right=True, top=True)

def plot_thermal_conductivity():
    """
    绘制热导率随时间变化的图。
    """
    # 加载热导率数据
    kappa = load_kappa()
    print(kappa.keys())  # 输出数据的键值信息

    # 计算运行平均值
    t = np.arange(1, kappa['kxi'].shape[0] + 1) * 0.001  # ns
    kappa['kyi_ra'] = running_ave(kappa['kyi'], t)
    kappa['kyo_ra'] = running_ave(kappa['kyo'], t)
    kappa['kxi_ra'] = running_ave(kappa['kxi'], t)
    kappa['kxo_ra'] = running_ave(kappa['kxo'], t)
    kappa['kz_ra'] = running_ave(kappa['kz'], t)

    # 创建绘图
    plt.figure(figsize=(12, 10))

    # 子图1: kyi
    plt.subplot(2, 2, 1)
    set_fig_properties([plt.gca()])
    plt.plot(t, kappa['kyi'], color='C7', alpha=0.5)
    plt.plot(t, kappa['kyi_ra'], linewidth=2)
    plt.xlim([0, 10])
    plt.gca().set_xticks(range(0, 11, 2))
    plt.ylim([-10, 10])
    plt.xlabel('time (ns)')
    plt.ylabel(r'$\kappa_{in}$ W/m/K')
    plt.title('(a)')

    # 子图2: kyo
    plt.subplot(2, 2, 2)
    set_fig_properties([plt.gca()])
    plt.plot(t, kappa['kyo'], color='C7', alpha=0.5)
    plt.plot(t, kappa['kyo_ra'], linewidth=2, color='C3')
    plt.xlim([0, 10])
    plt.gca().set_xticks(range(0, 11, 2))
    plt.ylim([-10, 10])
    plt.xlabel('time (ns)')
    plt.ylabel(r'$\kappa_{out}$ (W/m/K)')
    plt.title('(b)')

    # 子图3: kyi_ra + kyo_ra
    plt.subplot(2, 2, 3)
    set_fig_properties([plt.gca()])
    plt.plot(t, kappa['kyi_ra'], linewidth=2)
    plt.plot(t, kappa['kyo_ra'], linewidth=2, color='C3')
    plt.plot(t, kappa['kyi_ra'] + kappa['kyo_ra'], linewidth=2, color='k')
    plt.xlim([0, 10])
    plt.gca().set_xticks(range(0, 11, 2))
    plt.ylim([-10, 10])
    plt.xlabel('time (ns)')
    plt.ylabel(r'$\kappa$ (W/m/K)')
    plt.legend(['in', 'out', 'total'])
    plt.title('(c)')

    # 子图4: kyi_ra + kyo_ra, kxi_ra + kxo_ra, kz_ra
    plt.subplot(2, 2, 4)
    set_fig_properties([plt.gca()])
    plt.plot(t, kappa['kyi_ra'] + kappa['kyo_ra'], color='k', linewidth=2)
    plt.plot(t, kappa['kxi_ra'] + kappa['kxo_ra'], color='C0', linewidth=2)
    plt.plot(t, kappa['kz_ra'], color='C3', linewidth=2)
    plt.xlim([0, 10])
    plt.gca().set_xticks(range(0, 11, 2))
    plt.ylim([-10, 10])
    plt.xlabel('time (ns)')
    plt.ylabel(r'$\kappa$ (W/m/K)')
    plt.legend(['yy', 'xy', 'zy'])
    plt.title('(d)')

    plt.tight_layout()
    plt.savefig('./HNEMD.png', bbox_inches='tight')

def main():
    """
    主函数，调用绘图函数。
    """
    plot_thermal_conductivity()

if __name__ == '__main__':
    main()
