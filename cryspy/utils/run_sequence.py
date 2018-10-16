"""Defines the Ewa object which interfaces with Ewald"""
import os
import subprocess

import cryspy.io.edit_file as ef
import cryspy.io.read_file as rf

from cryspy.scripts.assign_charges import assign_charges

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
        self.output_file = open("prep.out", "w")
        if self.inputs["ewald"]:
            pref = "ew_"
        else:
            pref = "noew_"
        if self.inputs["self_consistent"]:
            post = "sc"
        else:
            post = "nosc"
        self.mode = pref + post
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
            high_level_pop_mol = rf.mol_from_gauss("high_pop_file", pop=self.inputs["high_pop_method"])
            shell_high.populate(high_level_pop_mol)
        else:
            shell_high = self.cell.make_cluster(self.inputs["clust_rad"])
            for atom in self.region_1:
                if atom in shell_high:
                    shell_high.remove(atom)
        low_level_pop_mol = rf.mol_from_gauss("low_pop_file", pop=self.inputs["low_pop_method"])
        shell_low = shell_high.copy()
        shell_low.populate(low_level_pop_mol)
        return shell_high, shell_low

    def run_ewald(self):
        if not os.path.exists(self.ewald_path):
            os.makedirs(self.ewald_path)
        os.chdir(self.ewald_path)
        # no stdout
        FNULL = open(os.devnull, 'w')

        ef.write_uc(self.inputs["name"] + ".uc", self.cell_vectors, self.inputs["aN"], self.inputs["bN"], self.inputs["cN"], self.region_1)
        ef.write_qc(self.inputs["name"] + ".qc", self.region_1)
        ef.write_ew_in(self.inputs["name"], "ewald.in." + self.inputs["name"], self.inputs["nchk"], self.inputs["nat"])
        ef.write_seed()
        # run Ewald
        subprocess.call("$FRO_EWALD < ewald.in." + self.inputs["name"], stdout=FNULL, shell=True)
        os.chdir(self.here)

        points = rf.read_points(self.inputs["name"]+"pts-fro")
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
        region_2, high_points = run_types[self.mode]
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
        self.self_consistent(None)
        ew_points = self.run_ewald()

        return region_2_low, ew_points

    def single_sc_loop(self, sc_loop, region_2):
        """Run a single iteration of the sc loop, with or without Ewald"""
        sc_name = "sc_" + self.inputs["name"]
        sc_loop += 1
        # Initial charges in mol
        old_charges = self.region_1.charges()

        if self.mode == "ewsc":
            points = self.run_ewald()
            ef.write_gauss(self.calc_name + ".com", self.region_1, points, self.calc_name + ".temp")
        else:
            ef.write_gauss(self.calc_name + ".com", self.region_1, region_2, self.calc_name + ".temp")

        subprocess.Popen("${FRO_GAUSS} " + self.calc_name + ".com", shell=True)
        # Calculate new charges

        intact_charges, new_energy, char_self, char_int = rf.read_g_char(sc_name + ".log", self.inputs["high_pop_method"], debug=True)

        # Correct charges if they are not perfectly neutral
        if sum(intact_charges) != 0.0:
            temp_correct = sum(intact_charges) / len(intact_charges)
            intact_charges = [i - temp_correct for i in intact_charges]

        # Damp the change in charges
        new_charges = [new * (1 - self.inputs["damping"]) + old * self.inputs["damping"] for new, old in zip(intact_charges, old_charges)]

        # Correct charges again (due to damping)
        if sum(new_charges) != 0.0:
            temp_correct = sum(new_charges) / len(new_charges)
            new_charges = [i - temp_correct for i in new_charges]
        # assign damped charges
        self.region_1.assign_charges(new_charges)

        if self.mode == "noew_sc":
            assign_charges(self.region_1, region_2)

        # Calculate deviation between initial and new charges
        deviation = sum([abs(i - j)
                         for (i, j) in zip(intact_charges, old_charges)]) / len(self.region_1)

        out_str = ("Iteration:", sc_loop, "Deviation:",
                   deviation, "Energy:", new_energy, "Charge self energy:", char_self, "Total - charge self:", new_energy - char_self)
        self.output_file.write(
            ("{:<6} {:<5} {:<6} {:10.6f} {:<6} {:10.6f} {:<6} {:10.6f} {:<6} {:10.6f}\n".format(*out_str)))
        self.output_file.flush()

        return deviation

    def self_consistent(self, region_2_high):
        """Run single iterations until the charge deviation is below the tol"""
        sc_iter = 0
        dev  = float("str")
        while dev > self.inputs["dev_tol"]:
            dev = self.single_sc_loop(sc_iter, region_2_high)
        return
