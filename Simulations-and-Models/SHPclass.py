"""
New, better optimised version
File containing classes for SH project mRNA translation
Running this does nothing on its own
"""

import random
import numpy as np
import matplotlib.pyplot as plt
import csv


class Site(object):
    # class for one site on lattice
    def __init__(self, k):
        # transition rate, w, for particular site, will be defined by codon in future
        self.k = k
        # whether site is occupied (0 or 1) by first point of particle not all 10, all initially unoppupied
        self.occ = 0


class Translation(object):
    # class for full lattice of sites
    def __init__(self, alpha, beta, rates_list, true_length, part_size, sweeps, gamma):
        # impur parameters
        self.true_length = true_length
        self.part_size = part_size
        self.sweeps = sweeps
        self.gamma = gamma
        self.alpha = alpha
        self.beta = beta

        # setting hopping rates
        ks_full_list = [alpha] + rates_list
        #ks_full_list.append(beta)
        ks_full_list[-1] = beta
        ks_full_list.append(0)

        self.lattice = np.empty(self.true_length + 2, dtype=object)
        for i in range(len(self.lattice)):
            self.lattice[i] = Site(ks_full_list[i])
        # set particle waiting to come onto lattice
        self.lattice[0].occ = 1

        # for calculating densities and current
        self.ns = np.zeros(len(self.lattice))
        self.qs = np.zeros(len(self.lattice))
        self.num_rib = 0

        # tracking time
        self.total_t = 0
        self.hit_ss = False
        self.t0 = 0
        self.loops = 0
        self.ss_frames = 0


    def get_occs_array(self):
        # get list of 1's and 0's for lattice including imaginary first and last sites
        occs_list = []
        for i in range(len(self.lattice)):
             occs_list.append(self.lattice[i].occ)
        return np.array(occs_list)


    def get_true_occs_array(self):
        # get list of 1's and 0's for lattice excluding imaginary first and last sites
        occs_list = []
        for i in range(1, len(self.lattice) - 1):
             occs_list.append(self.lattice[i].occ)
        return np.array(occs_list)


    def can_move(self, i):
        # check if particle on site i has space ahead to move into
        # scan ahead by particle size
        if i == 0:
            space = self.lattice[i + 1:i + 1 + self.part_size]
            occs = []
            for site in space:
                occs.append(site.occ)
            if 1 not in occs:
                move_possiblity = True
            else:
                move_possiblity = False

        elif i <= len(self.lattice) - 2 - self.part_size:
            if self.lattice[i + self.part_size].occ == 0:
                move_possiblity = True
            else:
                move_possiblity = False

        else:
            move_possiblity = True

        return move_possiblity


    def next_move(self):
        # go through one iteration by choosing next state out of different possibilities
        r = 0
        marker = 0
        new_state_probs = []
        # makes tuples with (site index of particle that can move, start, finish of where that slice of w lies on line 0 to r)
        for i in range(len(self.lattice) - 1):
            if self.lattice[i].occ == 1:
                if self.can_move(i):
                    new_state_probs.append((i, marker, marker + self.lattice[i].k, False))
                    marker += self.lattice[i].k
                    r += self.lattice[i].k
                else:
                    new_state_probs.append((i, marker, marker + self.gamma, True))
                    marker += self.gamma
                    r += self.gamma

        # where on line 0 to r to select next state
        state_choice_num = random.uniform(0, r)
        for state in new_state_probs:
            if state[1] <= state_choice_num <= state[2]:
                moving_site = state[0]
                collision = state[3]

        self.lattice[moving_site].occ = 0
        if collision:
            self.lattice[moving_site + self.part_size].occ = 0
        self.lattice[moving_site + 1].occ = 1

        # new particle ready to enter from left imaginary site
        self.lattice[0].occ = 1
        # particle off of lattice vanishes
        self.lattice[len(self.lattice) - 1].occ = 0

        # time for this is sampled from exp distribution
        delta_t = np.random.exponential(1/r)
        self.total_t += delta_t

        if self.hit_ss == True:
            self.ss_frames += 1
            self.ns += self.get_occs_array() * delta_t
            self.qs[moving_site] += 1
            self.num_rib += np.sum(self.get_true_occs_array())


    def check_ss(self, num_parts1, num_parts2):
        # compare mean number of particles from successive sets of steps and return whethr ss reached
        mean_num1 = sum(num_parts1) / len(num_parts1)
        mean_num2 = sum(num_parts2) / len(num_parts2)
        percent_diff = 100 * (abs(mean_num2 - mean_num1) / mean_num2)

        thresh = 2
        if percent_diff < thresh:
            self.hit_ss = True
            self.t0 = self.total_t


    def play(self):
        # play the simulation
        rows_written = 0
        segment = 1000
        num_parts1 = []
        num_parts2 = []
        with open('5e4_snaps_alpha60_beta60.csv', 'w', encoding='UTF8', newline='') as f:
            writer = csv.writer(f)
            oncemore = iter([True, False])
            while self.loops == 0 or not self.hit_ss or next(oncemore):
                for i in range(self.sweeps * self.true_length):
                    if not self.hit_ss:
                        if len(num_parts1) < segment:
                            num_parts1.append(self.get_true_occs_array().sum())
                        else:
                            num_parts2.append(self.get_true_occs_array().sum())
                        if len(num_parts1) == segment and len(num_parts2) == segment:
                            self.check_ss(num_parts1, num_parts2)
                            num_parts1 = []
                            num_parts2 = []
                    self.next_move()
                    if self.hit_ss == True and rows_written < 5*10**4:
                        random_0_to_1 = random.random()
                        if random_0_to_1 < 0.05:
                            writer.writerow(list(self.get_true_occs_array()))
                            rows_written += 1
                        if rows_written == 5*10**4:
                            print('csv done, rows written: ', rows_written)
                self.loops += 1


    def get_densities(self):
        # return densities at each site
        true_ns = self.ns[1:-1]
        densities = true_ns / (self.total_t - self.t0)

        return densities


    def get_currents(self):
        # return currents at each site and protein rate (current from final site)
        true_qs = self.qs[1:-1]
        currents = true_qs / (self.total_t - self.t0)
        protein_rate = currents[-1]

        return currents, protein_rate


    def make_plots(self, densities, currents):
        # densities and currents

        plt.figure(figsize=(12, 4))

        ax1 = plt.subplot2grid((3,7), (0,0), rowspan=3, colspan=3)
        ax2 = plt.subplot2grid((3,7), (0,4), rowspan=3, colspan=3)

        ax1.set_xlabel('i')
        ax2.set_xlabel('i')
        ax1.set_ylabel('Density, ' + r'$\rho_i$')
        ax2.set_ylabel('Current, J')

        # densities always between 0 and 1 so fix y axis limits
        ax1.set_ylim(0,1)

        ax1.plot(densities, c='blue')
        ax2.plot(currents, c='green')

        plt.suptitle(r'$\alpha$ = ' + str(self.alpha) + ', ' + r'$\beta$ = ' + str(self.beta))
        plt.show()
