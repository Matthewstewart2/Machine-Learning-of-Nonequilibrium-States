"""
File containing simulation for SH project mRNA translation
Produces a visual grid and other plots
"""

from SHPCorrclass import Translation
import numpy as np
import matplotlib.pyplot as plt
import csv

class Simulation(object):
    # class just runs simulation
    def run1(self, writer, sample_diff):
        # basic test of model
        alpha = 0.75
        beta = 0.75
        true_length = 100
        rates_list = [1] * true_length
        part_size = 1
        sweeps = 100
        gamma = 0

        lattice_model = Translation(alpha, beta, rates_list, true_length, part_size, sweeps, gamma)
        lattice_model.play(writer, sample_diff)
        #densities = lattice_model.get_densities()
        #currents, protein_rate = lattice_model.get_currents()
        #lattice_model.make_plots(densities, currents)


    def gamma_test(self):
        # simulation including real sequences and non-zero gamma
        rates_list = []
        filename = 'YDL184C_rates.DAT'
        with open(filename) as tsv:
            for line in csv.reader(tsv, delimiter=' '):
                rates_list.append(float(line[-1]))

        alpha_max = 10
        num_points = 10
        alphas = np.linspace(0.1, alpha_max, num_points)
        beta = rates_list[-1]
        true_length = len(rates_list)
        part_size = 9
        sweeps = 1000
        gammas = [0.1, 1.0, 10.0]
        plot_colours = ['.k', '.b', '.r']
        labels = [r'$\gamma$ = 0.1', r'$\gamma$ = 1.0', r'$\gamma$ = 10.0']

        for n in range(3):
            protein_alpha_rates = np.zeros(len(alphas))
            rates_err = np.zeros(len(protein_alpha_rates))
            for i in range(len(alphas)):
                protein_rate_results = np.zeros(10)
                for j in range(len(protein_rate_results)):
                    lattice_model = Translation(alphas[i], beta, rates_list, true_length, part_size, sweeps, gammas[n])
                    lattice_model.play()
                    protein_rate_results[j] = lattice_model.get_currents()[1]
                protein_alpha_rates[i] = np.mean(protein_rate_results)
                std_protein = np.std(protein_rate_results)
                rates_err[i] = std_protein / (len(protein_rate_results) ** 0.5)

            plt.errorbar(alphas, protein_alpha_rates, yerr=rates_err, fmt=plot_colours[n], capsize=3, label=labels[n])


        plt.xlabel(r'$\alpha$')
        plt.ylabel('Protein rate')
        plt.legend()
        plt.title('Protein Rate Against Initiation Rate')
        plt.show()


    def real_seq(self):
        # test of real sequences read in from file
        rates_list = []
        filename = 'YDL184C_rates.DAT'
        with open(filename) as tsv:
            for line in csv.reader(tsv, delimiter=' '):
                rates_list.append(float(line[-1]))

        alpha = 1.24654
        beta = rates_list[-1]
        true_length = len(rates_list)
        part_size = 9
        sweeps = 6000
        gamma = 0

        current_results = np.zeros(10)
        density_results = np.zeros(10)
        new_density_results = np.zeros(10)
        for i in range(len(current_results)):
            lattice_model = Translation(alpha, beta, rates_list, true_length, part_size, sweeps, gamma)
            lattice_model.play()
            densities = lattice_model.get_densities()
            currents, protein_rate = lattice_model.get_currents()
            new_density = lattice_model.get_new_density()
            density_results[i] = np.mean(densities)
            current_results[i] = np.mean(currents)
            new_density_results[i] = new_density

        mean_rho = np.mean(density_results)
        std_rho = np.std(density_results)
        sem_rho = std_rho / (len(density_results) ** 0.5)

        mean_J = np.mean(current_results)
        std_J = np.std(current_results)
        sem_J = std_J / (len(current_results) ** 0.5)

        print('mean rho = ' + str(mean_rho) + ' error = ' + str(sem_rho))
        print('mean J = ' + str(mean_J) + ' error = ' + str(sem_J))

def main():

    sample_diff = 0.05
    test = Simulation()
    with open('MCRealTimeCorrelationsDiffpoint05.csv', 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        for i in range(5000):
            test.run1(writer, sample_diff)
            print(i + 1)
    #test.gamma_test()
    #test.real_seq()

main()
