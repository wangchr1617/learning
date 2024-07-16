from pylab import *

##set figure properties
aw = 1.5
fs = 16
lw = 1.5
ms = 6
font = {'size'   : fs}
matplotlib.rc('font', **font)
matplotlib.rc('axes' , lw=aw)

def set_fig_properties(ax_list):
    tl = 6
    tw = 1.5
    tlm = 3
    
    for ax in ax_list:
        ax.tick_params(which='major', length=tl, width=tw)
        ax.tick_params(which='minor', length=tlm, width=tw)
        ax.tick_params(which='both', axis='both', direction='out', right=False, top=False)


loss = loadtxt('./loss.out')
loss[:,0] = np.arange(1, len(loss) + 1)*100
print("We have run %s steps!"%loss[-1, 0])
energy_train = loadtxt('./energy_train.out')
force_train = loadtxt('./force_train.out')
virial_train = loadtxt('./virial_train.out')
#energy_test = loadtxt('./energy_test.out')
#force_test = loadtxt('./force_test.out')
#virial_test = loadtxt('./virial_test.out')

figure(figsize=(12, 10))
subplot(2, 2, 1)
set_fig_properties([gca()])
loglog(loss[:, 0], loss[:, 1],  ls="-", lw=lw, c = "C8", label="Total")
loglog(loss[:, 0], loss[:, 2],  ls="-", lw=lw, c = "C0", label=r"$L_{1}$")
loglog(loss[:, 0], loss[:, 3],  ls="-", lw=lw, c = "C1", label=r"$L_{2}$")
loglog(loss[:, 0], loss[:, 4],  ls="-", lw=lw, c = "C2", label="E-train")
loglog(loss[:, 0], loss[:, 5],  ls="-", lw=lw, c = "C3", label="F-train")
loglog(loss[:, 0], loss[:, 6],  ls="-", lw=lw, c = "C4", label="V-train")
#loglog(loss[:, 0], loss[:, 7],  ls="-", lw=lw, c = "C5", label="E-test")
#loglog(loss[:, 0], loss[:, 8],  ls="-", lw=lw, c = "C6", label="F-test")
#loglog(loss[:, 0], loss[:, 9],  ls="-", lw=lw, c = "C7", label="V-test")
xlim(1e2,)
ylim(5e-4, 2e0)
xlabel('Generation')
ylabel('Loss')
legend(loc="upper right",  
        ncol = 3, 
        frameon = True,
        fontsize=12,
        labelspacing = 0,
        columnspacing = 0)
title("(a)")


subplot(2, 2, 2)
set_fig_properties([gca()])
plot(energy_train[:, 1], energy_train[:, 0], 'o', c="C2", ms = ms, alpha=0.5, label="Train")
#plot(energy_test[:, 1], energy_test[:, 0], 'o', c="C5", ms = ms, alpha=0.5, label="Test")
plot([np.min(energy_train)-0.1, np.max(energy_train)+0.1], [np.min(energy_train)-0.1, np.max(energy_train)+0.1], c = "grey", lw = 1)
xlim([np.min(energy_train)-0.1, np.max(energy_train)+0.1])
ylim([np.min(energy_train)-0.1, np.max(energy_train)+0.1])
xlabel('DFT energy (eV/atom)')
ylabel('NEP energy (eV/atom)')
legend(loc="upper left")
title("(b)")


subplot(2, 2, 3)
set_fig_properties([gca()])
plot(force_train[:, 3], force_train[:, 0], 'o', c="C3", ms = ms, alpha=0.5, label="Train")
plot(force_train[:, 4:6], force_train[:, 1:3], 'o', c="C3", ms = ms, alpha=0.5)
#plot(force_test[:, 3], force_test[:, 0], 'o', c="C6", ms = ms, alpha=0.5, label="Test")
#plot(force_test[:, 4:6], force_test[:, 1:3], 'o', c="C6", ms = ms, alpha=0.5)
plot([np.min(force_train)-1, np.max(force_train)+1], [np.min(force_train)-1, np.max(force_train)+1], c = "grey", lw = 1)
xlim([np.min(force_train)-1, np.max(force_train)+1])
ylim([np.min(force_train)-1, np.max(force_train)+1])
xlabel(r'DFT force (eV/$\rm{\AA}$)')
ylabel(r'NEP force (eV/$\rm{\AA}$)')
legend(loc="upper left")
title("(c)")


subplot(2, 2, 4)
set_fig_properties([gca()])
plot(virial_train[:, 6], virial_train[:, 0], 'o', c="C4", ms = ms, alpha=0.5, label="Train")
plot(virial_train[:, 7:12], virial_train[:, 1:6], 'o', c="C4", ms = ms, alpha=0.5)
#plot(virial_test[:, 6], virial_test[:, 0], 'o', c="C7", ms = ms, alpha=0.5, label="Test")
#plot(virial_test[:, 7:12], virial_test[:, 1:6], 'o', c="C7", ms = ms, alpha=0.5)
plot([np.min(virial_train)-1, np.max(virial_train)+1], [np.min(virial_train)-1, np.max(virial_train)+1], c = "grey", lw = 1)
xlim([np.min(virial_train)-1, np.max(virial_train)+1])
ylim([np.min(virial_train)-1, np.max(virial_train)+1])
xlabel('DFT virial (eV/atom)')
ylabel('NEP virial (eV/atom)')
legend(loc="upper left")
title("(d)")

subplots_adjust(wspace=0.35, hspace=0.35)
savefig("RMSE.png", bbox_inches='tight')
