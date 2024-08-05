# 文件名: plot_rmse_comparison.py
# 运行方法: python plot_rmse_comparison.py
# 功能描述: 从损失数据文件和训练数据文件中读取数据，并创建四个子图来比较损失、能量、力和应力的训练结果。

import matplotlib.pyplot as plt
import numpy as np

# 设置matplotlib全局参数
plt.rcParams['font.size'] = 16
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# 设置绘图参数
aw = 1.5
fs = 16
lw = 1.5
ms = 6

def set_fig_properties(ax_list):
    """
    设置图形属性。

    参数:
        ax_list (list): Axes对象列表。
    """
    tl = 6  # 大刻度长度
    tw = 1.5  # 大刻度宽度
    tlm = 3  # 小刻度长度
    
    for ax in ax_list:
        ax.tick_params(which='major', length=tl, width=tw)
        ax.tick_params(which='minor', length=tlm, width=tw)
        ax.tick_params(which='both', axis='both', direction='out', right=False, top=False)

def load_data():
    """
    加载数据文件。

    返回:
        tuple: 包含损失、能量、力和应力数据的元组。
    """
    loss = np.loadtxt('./loss.out')
    loss[:, 0] = np.arange(1, len(loss) + 1) * 100
    print(f"We have run {loss[-1, 0]:.0f} steps!")

    energy_train = np.loadtxt('./energy_train.out')
    force_train = np.loadtxt('./force_train.out')
    virial_train = np.loadtxt('./virial_train.out')

    return loss, energy_train, force_train, virial_train

def plot_training_results():
    """
    绘制损失、能量、力和应力的比较图。
    """
    # 加载数据
    loss, energy_train, force_train, virial_train = load_data()

    plt.figure(figsize=(12, 10))

    # 绘制损失曲线
    plt.subplot(2, 2, 1)
    set_fig_properties([plt.gca()])
    plt.loglog(loss[:, 0], loss[:, 1], ls="-", lw=lw, c="C8", label="Total")
    plt.loglog(loss[:, 0], loss[:, 2], ls="-", lw=lw, c="C0", label=r"$L_{1}$")
    plt.loglog(loss[:, 0], loss[:, 3], ls="-", lw=lw, c="C1", label=r"$L_{2}$")
    plt.loglog(loss[:, 0], loss[:, 4], ls="-", lw=lw, c="C2", label="E-train")
    plt.loglog(loss[:, 0], loss[:, 5], ls="-", lw=lw, c="C3", label="F-train")
    plt.loglog(loss[:, 0], loss[:, 6], ls="-", lw=lw, c="C4", label="V-train")
    plt.xlim(1e2,)
    plt.ylim(5e-4, 2e0)
    plt.xlabel('Generation')
    plt.ylabel('Loss')
    plt.legend(loc="upper right", ncol=3, frameon=True, fontsize=12, labelspacing=0, columnspacing=0)
    plt.title("(a)")

    # 绘制能量散点图
    plt.subplot(2, 2, 2)
    set_fig_properties([plt.gca()])
    plt.plot(energy_train[:, 1], energy_train[:, 0], 'o', c="C2", ms=ms, alpha=0.5, label="Train")
    plt.plot([np.min(energy_train) - 0.1, np.max(energy_train) + 0.1], [np.min(energy_train) - 0.1, np.max(energy_train) + 0.1], c="grey", lw=1)
    plt.xlim([np.min(energy_train) - 0.1, np.max(energy_train) + 0.1])
    plt.ylim([np.min(energy_train) - 0.1, np.max(energy_train) + 0.1])
    plt.xlabel('DFT energy (eV/atom)')
    plt.ylabel('NEP energy (eV/atom)')
    plt.legend(loc="upper left")
    plt.title("(b)")

    # 绘制力散点图
    plt.subplot(2, 2, 3)
    set_fig_properties([plt.gca()])
    plt.plot(force_train[:, 3], force_train[:, 0], 'o', c="C3", ms=ms, alpha=0.5, label="Train")
    plt.plot(force_train[:, 4:6], force_train[:, 1:3], 'o', c="C3", ms=ms, alpha=0.5)
    plt.plot([np.min(force_train) - 1, np.max(force_train) + 1], [np.min(force_train) - 1, np.max(force_train) + 1], c="grey", lw=1)
    plt.xlim([np.min(force_train) - 1, np.max(force_train) + 1])
    plt.ylim([np.min(force_train) - 1, np.max(force_train) + 1])
    plt.xlabel(r'DFT force (eV/$\rm{\AA}$)')
    plt.ylabel(r'NEP force (eV/$\rm{\AA}$)')
    plt.legend(loc="upper left")
    plt.title("(c)")

    # 绘制应力散点图
    plt.subplot(2, 2, 4)
    set_fig_properties([plt.gca()])
    plt.plot(virial_train[:, 6], virial_train[:, 0], 'o', c="C4", ms=ms, alpha=0.5, label="Train")
    plt.plot(virial_train[:, 7:12], virial_train[:, 1:6], 'o', c="C4", ms=ms, alpha=0.5)
    plt.plot([np.min(virial_train) - 1, np.max(virial_train) + 1], [np.min(virial_train) - 1, np.max(virial_train) + 1], c="grey", lw=1)
    plt.xlim([np.min(virial_train) - 1, np.max(virial_train) + 1])
    plt.ylim([np.min(virial_train) - 1, np.max(virial_train) + 1])
    plt.xlabel('DFT virial (eV/atom)')
    plt.ylabel('NEP virial (eV/atom)')
    plt.legend(loc="upper left")
    plt.title("(d)")

    plt.subplots_adjust(wspace=0.35, hspace=0.35)
    plt.savefig("RMSE.png", bbox_inches='tight')

def main():
    """
    主函数，调用绘图函数。
    """
    plot_training_results()

if __name__ == '__main__':
    main()
