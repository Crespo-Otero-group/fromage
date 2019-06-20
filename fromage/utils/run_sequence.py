"""Defines the Ewa object which interfaces with Ewald"""
import os
import subprocess
import time
import sys

import fromage.io.edit_file as ef
import fromage.io.read_file as rf

from fromage.scripts.fro_assign_charges import assign_charges

class RunSeq(object):
    """
    Class which sets up the order of operations for preparing calculation

    Attributes
    ----------
    region_1 : Mol object
        The atoms in the central molecule
    cell : Mol object
        The atoms in the unit cell
    inputs : dict
        The input keywords
    mode : str
        String summarising the kind of run sequence required. options are:
        noew_nosc : EC
        noew_sc : SC-EC
        ew_nosc : EEC
        ew_sc : SC-EEC
    """
    def __init__(self, region_1, cell, inputs):
        self.region_1 = region_1
        self.cell = cell
        self.inputs = inputs
        if self.inputs["ewald"]:
            pref = "ew_"
        else:
            pref = "noew_"
        if self.inputs["self_consistent"]:
            post = "sc"
        else:
            post = "nosc"
        self.mode = pref + post
        # dirs
        self.here = os.getcwd()
        self.ewald_path = os.path.join(self.here,"ewald/")
        self.out_file = open("prep.out","a")
        return

    def write_out(self,string):
        self.out_file.write(string)
        self.out_file.flush()
        return

    def make_region_2(self):
        """
        Get region 2 Mols with different charges

        Returns
        -------
        shell_high : Mol object
            Region 2 molecules with high level of theory charges
        shell_low : Mole object
            Region 2 molecules with low level of theory charges

        """
        if self.inputs["target_shell"]:
            shell_high = rf.mol_from_file(self.inputs["target_shell"])
            self.write_out("Outer region read in with " + str(len(shell_high)) + " atoms.\n")
            high_level_pop_mol = rf.mol_from_gauss(self.inputs["high_pop_file"], pop=self.inputs["high_pop_method"])
            shell_high.populate(high_level_pop_mol)
        else:
            shell_high = self.cell.make_cluster(self.inputs["clust_rad"], central_mol = self.region_1, mode = self.inputs["clust_mode"])
            for atom_i in self.region_1:
                for atom_j in shell_high:
                    if atom_i.very_close(atom_j):
                        shell_high.remove(atom_j)
                        break
            self.write_out("Outer region generated with " + str(len(shell_high)) + " atoms.\n")
        low_level_pop_mol = rf.mol_from_gauss(self.inputs["low_pop_file"], pop=self.inputs["low_pop_method"])
        shell_low = shell_high.copy()
        shell_low.populate(low_level_pop_mol)
        return shell_low, shell_high

    def run_ewald(self, calc_name=None):
        if calc_name == None:
            calc_name = self.inputs["name"]
        if not os.path.exists(self.ewald_path):
            os.makedirs(self.ewald_path)
        os.chdir(self.ewald_path)
        # no stdout
        FNULL = open(os.devnull, 'w')

        ef.write_uc(calc_name + ".uc", self.inputs["vectors"], self.inputs["an"], self.inputs["bn"], self.inputs["cn"], self.cell)
        ef.write_qc(calc_name + ".qc", self.region_1)
        ef.write_ew_in(calc_name, "ewald.in." + calc_name, self.inputs["nchk"], self.inputs["nat"])
        ef.write_seed()
        # run Ewald
        self.write_out("Ewald calculation started\n")
        ew_start = time.time()
        subprocess.call("${FRO_EWALD} < ewald.in." + calc_name, stdout=FNULL, shell=True)
        ew_end = time.time()
        self.write_out("Ewald calculation finished after "+str(round(ew_end - ew_start,3))+" s\n")
        points = rf.read_points(calc_name + ".pts-fro")
        if len(points) == 0:
            self.write_out("Something went wrong with the Ewald calculation, stopping...\n")
            sys.exit()
        os.chdir(self.here)

        return points

    def run(self):
        """
        Run the calculation for the corresponding self.mode

        Returns
        -------
        region_2 : Mol object
            Region 2 atoms with low level of theory charges
        high_points : Mol object
            Points that will embed mh, regardless of self.mode
        """
        run_types = {"noew_nosc":self.run_ec,
                    "noew_sc":self.run_scec,
                    "ew_nosc":self.run_eec,
                    "ew_sc":self.run_sceec}
        # execute the appropriate run type
        region_2, high_points = run_types[self.mode]()
        self.out_file.close()
        return region_2, high_points

    def run_ec(self):
        region_2_low , region_2_high = self.make_region_2()

        return region_2_low, region_2_high

    def run_scec(self):
        region_2_low , region_2_high = self.make_region_2()
        self.self_consistent(region_2_high)

        return region_2_low, region_2_high

    def run_eec(self):
        region_2_low , region_2_high = self.make_region_2()
        ew_points = self.run_ewald()

        return region_2_low, ew_points

    def run_sceec(self):
        region_2_low , region_2_high = self.make_region_2()
        self.self_consistent(None) # here, the None argument means that the initial background has yet to be computed
        ew_points = self.run_ewald()

        return region_2_low, ew_points

    def single_sc_loop(self, sc_loop, initial_bg):
        """Run a single iteration of the sc loop, with or without Ewald"""
        sc_name = "sc_" + self.inputs["name"]
        # Initial charges in mol
        old_charges = self.region_1.charges()

        # if sc_eec then there is no initial_bg so it needs to be computed
        if self.mode == "ew_sc":
            points = self.run_ewald(calc_name = sc_name)
            initial_bg = points

        ef.write_gauss(sc_name + ".com", self.region_1, initial_bg, self.inputs["sc_temp"])

        subprocess.call("${FRO_GAUSS} " + sc_name + ".com", shell=True)
        # Calculate new charges

        intact_charges, new_energy, char_self, char_int = rf.read_g_char(sc_name + ".log", self.inputs["high_pop_method"], debug=True)

        # Correct charges if they are not perfectly neutral
        if sum(intact_charges) != 0.0:
            temp_correct = sum(intact_charges) / len(intact_charges)
            intact_charges = [i - temp_correct for i in intact_charges]

        dummy_mol = self.region_1.copy()
        dummy_mol.raw_assign_charges(intact_charges)
        self.region_1.populate(dummy_mol)

        # Damp the change in charges
        new_charges = [new * (1 - self.inputs["damping"]) + old * self.inputs["damping"] for new, old in zip(self.region_1.charges(), old_charges)]

        # Correct charges again (due to damping)
        if sum(new_charges) != 0.0:
            temp_correct = sum(new_charges) / len(new_charges)
            new_charges = [i - temp_correct for i in new_charges]

        # assign damped charges
        self.region_1.raw_assign_charges(new_charges)
        self.cell.populate(self.region_1)

        if self.mode == "noew_sc":
            assign_charges(self.region_1, initial_bg)

        # Calculate deviation between initial and new charges
        deviation = sum([abs(i - j)
                         for (i, j) in zip(self.region_1.charges(), old_charges)]) / len(self.region_1)

        out_str = ("Iteration:", sc_loop, "Deviation:",
                   deviation, "Energy:", new_energy, "Charge self energy:", char_self, "Total - charge self:", new_energy - char_self)
        self.write_out("{:<6} {:<5} {:<6} {:10.6f} {:<6} {:10.6f} {:<6} {:10.6f} {:<6} {:10.6f}\n".format(*out_str))

        return deviation

    def self_consistent(self, initial_bg):
        """Run single iterations until the charge deviation is below the tol"""
        sc_iter = 0
        dev  = float("inf")
        while dev > self.inputs["dev_tol"]:
            sc_iter += 1
            dev = self.single_sc_loop(sc_iter, initial_bg)
        self.write_out("Tolerance reached: " + str(dev) + " < " + str(self.inputs["dev_tol"]) + "\n")
        return
