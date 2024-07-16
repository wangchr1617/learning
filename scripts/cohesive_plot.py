import matplotlib.pyplot as plt

data_dict = {}
with open('ev.txt', 'r') as file:
    label = None 
    for line in file:
        try:
            vol, ene = map(float, line.split(','))
            if label is not None:
                data_dict[label]['volumes'].append(vol)
                data_dict[label]['energies'].append(ene)
        except ValueError:
            label = line.strip()  
            data_dict[label] = {'volumes': [], 'energies': []} 

plt.figure(figsize=(10, 6), dpi=300)
for label, data in data_dict.items():
    plt.plot(data['volumes'], data['energies'], linestyle='--', alpha=0.7, label=label)
#    plt.scatter(data['volumes'], data['energies'], alpha=0.7)

plt.xlabel(r'Volume / $Å^3$')
plt.ylabel('Energy per atom/eV')
plt.legend(loc="best")
plt.savefig('./ev.png', bbox_inches='tight')

