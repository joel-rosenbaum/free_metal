import math
import json
import re

class IonAdjuster:
    def __init__(
        self,
        ph_value,
        ionic_strength,
        temperature_celsius,
        standard_ionic_strength=0.10,
        standard_temperature_celsius=20.0
    ):
        self.ph_value = ph_value
        self.ionic_strength = ionic_strength
        self.temperature_celsius = temperature_celsius
        self.standard_ionic_strength = standard_ionic_strength
        self.standard_temperature_celsius = standard_temperature_celsius
        self.standard_temperature_kelvin = self.standard_temperature_celsius + 273.15

    def compute_hydrogen_correction(self):
        """Compute activity-corrected [H+] based on pH and ionic strength."""
        debye_huckel_b = 0.522932 * math.exp(0.0327016 * self.temperature_celsius) + 4.015942
        gamma_h = (
            0.145045 * math.exp(-debye_huckel_b * self.ionic_strength)
            + 0.063546 * math.exp(-43.97704 * self.ionic_strength)
            + 0.695634
        )
        log_hydrogen = math.log10((10 ** (-self.ph_value)) / gamma_h)
        hydrogen_molar = 10 ** log_hydrogen
        return debye_huckel_b, gamma_h, log_hydrogen, hydrogen_molar

    def compute_dielectric_and_log_f(self, temperature_celsius, ionic_strength):
        """Compute dielectric constant and log_f for activity corrections."""
        temperature_kelvin = temperature_celsius + 273.15
        dielectric_constant = (
            87.7251
            - 0.3974762 * temperature_celsius
            + 0.0008253 * (temperature_celsius ** 2)
        )
        activity_coefficient_a = 1824600.0 / ((dielectric_constant * temperature_kelvin) ** 1.5)
        log_f_value = activity_coefficient_a * (
            (ionic_strength ** 0.5) / (1.0 + ionic_strength ** 0.5)
            - 0.25 * ionic_strength
        )
        return dielectric_constant, activity_coefficient_a, log_f_value

    def apply_ionic_strength_correction(self, log_k_standard, charge_factor_g, log_f_standard, log_f_experimental):
        """ logK_corrected = logK_standard + (charge_factor) * (log_f_standard - log_f_experimental)"""
        return log_k_standard + charge_factor_g * (log_f_standard - log_f_experimental)

    def apply_temperature_correction(self, log_k_ionic, enthalpy_kcal_per_mol, temperature_kelvin_std, temperature_kelvin_exp):
        """Apply approximate van't Hoff enthalpy correction."""
        inverse_factor = 1.0 / 0.0045765
        return log_k_ionic + enthalpy_kcal_per_mol * inverse_factor * (
            (1.0 / temperature_kelvin_std) - (1.0 / temperature_kelvin_exp)
        )

    def adjust_log_k_for_ionic_strength_and_temperature(
        self,
        log_k_standard,
        charge_factor_g,
        enthalpy_kcal_per_mol,
        log_f_standard,
        log_f_experimental,
        temperature_kelvin_std,
        temperature_kelvin_exp
    ):
        ionic_corrected = self.apply_ionic_strength_correction(
            log_k_standard, charge_factor_g, log_f_standard, log_f_experimental
        )
        temperature_corrected = self.apply_temperature_correction(
            ionic_corrected, enthalpy_kcal_per_mol, temperature_kelvin_std, temperature_kelvin_exp
        )
        return temperature_corrected
    # end of IonAdjuster


class ConstantUtils:
    def __init__(self, proton_json_path, chelator_json_path):
        with open(proton_json_path, "r") as f:
            self.protonation_data = json.load(f)
        with open(chelator_json_path, "r") as f:
            self.chelator_metal_data = json.load(f)

    def get_protonation_constants(self, chelator_name):
        block = None
        for entry in self.protonation_data["HYDROGEN_CONSTANTS"]:
            if entry["chelator"] == chelator_name:
                block = entry
                break
        if not block:
            raise ValueError(f"{chelator_name} not found in HYDROGEN_CONSTANTS")

        proton_list = []
        for i in range(1, 5):
            val_logk = block[f"H{i}"]
            val_dh = block[f"dH{i}"]
            if val_logk > 0.0:
                proton_list.append((val_logk, val_dh))
        return proton_list

    def get_metal_constants(self, chelator_name, metal_name):
        block = None
        for item in self.chelator_metal_data["CHELATOR_CONSTANTS"]:
            if item["chelator"] == chelator_name:
                block = item
                break

        if not block:
            raise ValueError(f"{chelator_name} not found in CHELATOR_CONSTANTS")

        if metal_name not in block["metals"]:
            raise ValueError(f"{metal_name} not found under {chelator_name}")

        metal_dict = block["metals"][metal_name]

        metal_ligand_log_k = metal_dict["ML_logK"]
        metal_ligand_enthalpy = metal_dict["ML_dH"]
        metal_h_ligand_log_k = metal_dict["MHL_logK"]
        metal_h_ligand_enthalpy = metal_dict["MHL_dH"]

        stoichiometry_string = metal_dict.get("stoic", "1:1")  # default stoichiometry if missing
        metal_stoich_str, ligand_stoich_str = stoichiometry_string.split(":")
        metal_stoichiometry = int(metal_stoich_str)
        ligand_stoichiometry = int(ligand_stoich_str)

        # Try to parse valence from something like "Ca2", "Mg2", etc.
        match = re.search(r'\d+', metal_name)
        if match:
            metal_valence = int(match.group())
        else:
            raise ValueError(f"Unable to extract valence from metal name '{metal_name}'")

        return {
            "metal_ligand_log_k": metal_ligand_log_k,
            "metal_ligand_enthalpy": metal_ligand_enthalpy,
            "metal_h_ligand_log_k": metal_h_ligand_log_k,
            "metal_h_ligand_enthalpy": metal_h_ligand_enthalpy,
            "metal_stoichiometry": metal_stoichiometry,
            "ligand_stoichiometry": ligand_stoichiometry,
            "metal_valence": metal_valence
        }
    # end of ConstantUtils


class EquilibriumUtils:
    @staticmethod
    def calc_charge_factor_for_protonation(num_sites):
        return 2.0 * num_sites

    @staticmethod
    def calc_charge_factor_for_metal_binding(chelator_valence, metal_valence):
        return 2.0 * chelator_valence * metal_valence

    @staticmethod
    def compute_apparent_binding_constant(log_k_adjusted, adjusted_proton_log_k_list, hydrogen_molar):
        """Compute pH-adjusted K from log10(K) plus protonation steps."""
        binding_k = 10 ** log_k_adjusted
        proton_k_values = [10 ** val for val in adjusted_proton_log_k_list]
        denominator = 1.0
        product_accumulator = 1.0
        for proton_k in proton_k_values:
            product_accumulator *= (proton_k * hydrogen_molar)
            denominator += product_accumulator
        return binding_k / denominator
    
    @staticmethod
    def solve_free_metal(
        metals_total,
        chelators_total,
        adjusted_logk,      # 2D matrix [m_idx][c_idx]
        metal_stoich,       # 2D matrix [m_idx][c_idx]
        ligand_stoich,      # 2D matrix [m_idx][c_idx]
        max_iter=500,
        tolerance=1e-10
    ):
        """
        Iteratively solve for free metal and chelator concentrations.

        Returns:
            (metals_free, chelators_free, iterations)
        """
        num_metals = len(metals_total)
        num_chelators = len(chelators_total)

        metals_free = [mtot * 0.5 for mtot in metals_total]
        chelators_free = [ctot * 0.5 for ctot in chelators_total]

        iterations = 0

        for _ in range(max_iter):
            iterations += 1
            converged = True  

            # Update chelator
            for c_idx in range(num_chelators):
                c_tot = chelators_total[c_idx]
                if c_tot <= 0:
                    chelators_free[c_idx] = 0.0
                    continue

                sum_binding = 0.0
                for m_idx in range(num_metals):
                    x = metal_stoich[m_idx][c_idx]
                    y = ligand_stoich[m_idx][c_idx]
                    if (x == 0 and y == 0):
                        continue

                    logk_val = adjusted_logk[m_idx][c_idx]
                    K_val = 10 ** logk_val
                    sum_binding += y * K_val * (metals_free[m_idx] ** x) * (chelators_free[c_idx] ** (y - 1))

                new_c_free = c_tot / (1.0 + sum_binding)
                diff = abs(new_c_free - chelators_free[c_idx])
                if diff >= tolerance * max(1e-9, chelators_free[c_idx]):
                    converged = False
                chelators_free[c_idx] = new_c_free

            # Update metal
            for m_idx in range(num_metals):
                m_tot = metals_total[m_idx]
                if m_tot <= 0:
                    metals_free[m_idx] = 0.0
                    continue

                sum_binding = 0.0
                for c_idx in range(num_chelators):
                    x = metal_stoich[m_idx][c_idx]
                    y = ligand_stoich[m_idx][c_idx]
                    if (x == 0 and y == 0):
                        continue

                    logk_val = adjusted_logk[m_idx][c_idx]
                    K_val = 10 ** logk_val
                    sum_binding += x * K_val * (metals_free[m_idx] ** (x - 1)) * (chelators_free[c_idx] ** y)

                new_m_free = m_tot / (1.0 + sum_binding)
                diff = abs(new_m_free - metals_free[m_idx])
                if diff >= tolerance * max(1e-9, metals_free[m_idx]):
                    converged = False
                metals_free[m_idx] = new_m_free

            if converged:
                break

        return metals_free, chelators_free, iterations

    @staticmethod
    def adjust_protonation_constants(
        ion_adjuster,
        protonation_steps,     # list of (logK, dH)
        log_f_std, 
        log_f_exp,
        temp_kelvin_std,
        temp_kelvin_exp
    ):
        """Adjust each proton site logK for ionic strength & temperature."""
        adjusted_list = []
        total_proton_sites = len(protonation_steps)
        for idx, (proto_logk, proto_dh) in enumerate(protonation_steps):
            sites_remaining = total_proton_sites - idx
            g_factor_proton = EquilibriumUtils.calc_charge_factor_for_protonation(sites_remaining)

            corrected = ion_adjuster.adjust_log_k_for_ionic_strength_and_temperature(
                proto_logk,
                g_factor_proton,
                proto_dh,
                log_f_std,
                log_f_exp,
                temp_kelvin_std,
                temp_kelvin_exp
            )
            adjusted_list.append(corrected)
        return adjusted_list

    @staticmethod
    def adjust_metal_constants(
        ion_adjuster,
        ml_logk,           # ML_logK (standard)
        ml_dh,             # ML enthalpy
        mlh_logk,          # MHL_logK (standard)
        mlh_dh,            # MHL enthalpy
        metal_valence,
        chelator_valence,
        log_f_std,
        log_f_exp,
        temp_kelvin_std,
        temp_kelvin_exp
    ):
        """Adjust ML and MHL standard logK for ionic strength & temperature."""
        g_factor_metal = EquilibriumUtils.calc_charge_factor_for_metal_binding(
            chelator_valence, metal_valence
        )

        corr_ml = ion_adjuster.adjust_log_k_for_ionic_strength_and_temperature(
            ml_logk,
            g_factor_metal,
            ml_dh,
            log_f_std,
            log_f_exp,
            temp_kelvin_std,
            temp_kelvin_exp
        )

        # If metal-hydrogen-ligand logK is zero or missing, just skip
        if mlh_logk == 0:
            corr_mlh = 0.0
        else:
            corr_mlh = ion_adjuster.adjust_log_k_for_ionic_strength_and_temperature(
                mlh_logk,
                g_factor_metal,
                mlh_dh,
                log_f_std,
                log_f_exp,
                temp_kelvin_std,
                temp_kelvin_exp
            )

        return corr_ml, corr_mlh
        # end of EquilibriumUtils
        
class Calculator:
    @staticmethod
    def mainblock(
        const_utils,
        chelators,
        metals,
        ph_value,
        ionic_strength_value,
        temp_celsius
    ):
        """Static method to run the entire calculation workflow and return results as dictionary."""
        ion_adj = IonAdjuster(
            ph_value,
            ionic_strength_value,
            temp_celsius
        )

        # Compute hydrogen correction
        _, _, _, hydrogen_molar = ion_adj.compute_hydrogen_correction()

        # Compute log_f values at standard & experimental conditions
        _, _, log_f_std_val = ion_adj.compute_dielectric_and_log_f(
            ion_adj.standard_temperature_celsius, 
            ion_adj.standard_ionic_strength
        )
        _, _, log_f_exp_val = ion_adj.compute_dielectric_and_log_f(
            temp_celsius, 
            ionic_strength_value
        )

        temp_kelvin_std = ion_adj.standard_temperature_kelvin
        temp_kelvin_exp = temp_celsius + 273.15

        # Initialize data structures
        num_metals = len(metals)
        num_chelators = len(chelators)

        adjusted_logk_matrix = [[] for _ in range(num_metals)]
        metal_stoich_matrix  = [[] for _ in range(num_metals)]
        ligand_stoich_matrix = [[] for _ in range(num_metals)]

        # Loop over each metal–chelator pair
        for m_idx, metal in enumerate(metals):
            metal_name = metal["name"]
            for c_idx, chelator in enumerate(chelators):
                chelator_name = chelator["name"]

                # Standard constants
                cdata = const_utils.get_metal_constants(chelator_name, metal_name)
                ml_logk  = cdata["metal_ligand_log_k"]
                ml_dh    = cdata["metal_ligand_enthalpy"]
                mlh_logk = cdata["metal_h_ligand_log_k"]
                mlh_dh   = cdata["metal_h_ligand_enthalpy"]
                metal_valence = cdata["metal_valence"]

                protonation_steps = const_utils.get_protonation_constants(chelator_name)
                chelator_valence = len(protonation_steps)
                
                # Adjust metal binding for ion activity
                corr_ml, corr_mlh = EquilibriumUtils.adjust_metal_constants(
                    ion_adjuster = ion_adj,
                    ml_logk      = ml_logk,
                    ml_dh        = ml_dh,
                    mlh_logk     = mlh_logk,
                    mlh_dh       = mlh_dh,
                    metal_valence    = metal_valence,
                    chelator_valence = chelator_valence,
                    log_f_std    = log_f_std_val,
                    log_f_exp    = log_f_exp_val,
                    temp_kelvin_std = temp_kelvin_std,
                    temp_kelvin_exp = temp_kelvin_exp
                )

                # Adjust protonation for ion activity
                adjusted_proton_log_k_list = EquilibriumUtils.adjust_protonation_constants(
                    ion_adjuster    = ion_adj,
                    protonation_steps = protonation_steps,
                    log_f_std       = log_f_std_val,
                    log_f_exp       = log_f_exp_val,
                    temp_kelvin_std = temp_kelvin_std,
                    temp_kelvin_exp = temp_kelvin_exp
                )

                # Compute apparent logK with all factors included
                apparent_ml = EquilibriumUtils.compute_apparent_binding_constant(
                    corr_ml, adjusted_proton_log_k_list, hydrogen_molar
                )

                if mlh_logk != 0:
                    apparent_mlh = EquilibriumUtils.compute_apparent_binding_constant(
                        corr_mlh, adjusted_proton_log_k_list[1:], hydrogen_molar
                    )
                    net_logk = math.log10(apparent_ml + apparent_mlh)
                else:
                    net_logk = math.log10(apparent_ml)

                adjusted_logk_matrix[m_idx].append(net_logk)
                metal_stoich_matrix[m_idx].append(cdata["metal_stoichiometry"])
                ligand_stoich_matrix[m_idx].append(cdata["ligand_stoichiometry"])

        # Solve equilibrium
        metals_total    = [m["concentration"] for m in metals]
        chelators_total = [c["concentration"] for c in chelators]

        metals_free, chelators_free, iters = EquilibriumUtils.solve_free_metal(
            metals_total,
            chelators_total,
            adjusted_logk_matrix,
            metal_stoich_matrix,
            ligand_stoich_matrix
        )

        # Compute final Kd for each pair
        kd_values = []
        for m_idx, metal in enumerate(metals):
            for c_idx, chelator in enumerate(chelators):
                net_logk = adjusted_logk_matrix[m_idx][c_idx]
                kd = 1.0 / (10 ** net_logk)
                kd_values.append({
                    "metal": metal["name"],
                    "chelator": chelator["name"],
                    "Kd": kd
                })

        # Return all numerical results in a dictionary
        return {
            "metals_free": metals_free,
            "chelators_free": chelators_free,
            "kd_values": kd_values,
            "iterations": iters
        }
    

if __name__ == "__main__":
    const_utils = ConstantUtils("webmaxc_protonation.json", "webmaxc_chelator.json")

    chelators = [
        {"name": "EGTA", "concentration": .005},
    ]
    metals = [
        {"name": "Ca2", "concentration": .006},
        {"name": "Zn2", "concentration": .001}
    ]

    ph_value = 8.2
    ionic_strength_value = 0.15
    temp_celsius = 23.0

    results = Calculator.mainblock(
        const_utils   = const_utils,
        chelators     = chelators,
        metals        = metals,
        ph_value      = ph_value,
        ionic_strength_value = ionic_strength_value,
        temp_celsius  = temp_celsius
    )

    metals_free = results["metals_free"]
    chelators_free = results["chelators_free"]
    kd_values = results["kd_values"]

    print(f"Converged in {results['iterations']} iterations.\n")

    for idx, metal in enumerate(metals):
        m_free = metals_free[idx]
        print(f"Free metal {metal['name']} = {m_free:.6e} M")

    for idx, chelator in enumerate(chelators):
        c_free = chelators_free[idx]
        print(f"Free chelator {chelator['name']} = {c_free:.6e} M")

    print("")
    for kv in kd_values:
        print(f"Kd for {kv['chelator']}-{kv['metal']}: {kv['Kd']:.6e} M")