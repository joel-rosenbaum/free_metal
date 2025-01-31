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
        return log_k_standard + charge_factor_g * (log_f_standard - log_f_experimental)

    def apply_temperature_correction(self, log_k_ionic, enthalpy_kcal_per_mol, temperature_kelvin_std, temperature_kelvin_exp):
        inverse_factor = 1.0 / 0.0045765
        return log_k_ionic + enthalpy_kcal_per_mol * inverse_factor * ((1.0 / temperature_kelvin_std) - (1.0 / temperature_kelvin_exp))

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

        stoichiometry_string = metal_dict.get("stoic", "1:1") # default stoichiometry is 1:1 if missing
        metal_stoich_str, ligand_stoich_str = stoichiometry_string.split(":")
        metal_stoichiometry = int(metal_stoich_str)
        ligand_stoichiometry = int(ligand_stoich_str)
        
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
        """
        Combines chelator's protonation steps with a single metal–ligand logK 
        to produce an 'apparent' binding constant at the given [H+].
        """
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
        adjusted_logk,      # Single net log K (shape [m_idx][c_idx])
        metal_stoich,       # M^x
        ligand_stoich,      # L^y
        max_iter=500,
        tolerance=1e-10
    ):
        """
        Iteratively solve for free metal and chelator concentrations.

        Returns:
            metals_free: List of free metal concentrations.
            chelators_free: List of free chelator concentrations.
            iterations: Number of iterations taken to converge.
        """
        num_metals = len(metals_total)
        num_chels = len(chelators_total)

        metals_free = [mtot * 0.5 for mtot in metals_total]
        chelators_free = [ctot * 0.5 for ctot in chelators_total]

        iterations = 0  # Initialize iteration counter

        for iteration in range(max_iter):
            iterations += 1  # Increment iteration counter
            converged = True  

            # Update chelator concentration
            for c_idx in range(num_chels):
                c_tot = chelators_total[c_idx]
                if c_tot <= 0:
                    chelators_free[c_idx] = 0.0
                    continue

                sum_binding = 0.0
                for m_idx in range(num_metals):
                    x = metal_stoich[m_idx][c_idx]
                    y = ligand_stoich[m_idx][c_idx]
                    if x == 0 and y == 0:
                        continue

                    logk_val = adjusted_logk[m_idx][c_idx]
                    K_val = 10 ** logk_val
                    sum_binding += y * K_val * (metals_free[m_idx] ** x) * (chelators_free[c_idx] ** (y - 1))

                new_c_free = c_tot / (1.0 + sum_binding)
                diff = abs(new_c_free - chelators_free[c_idx])
                if diff < tolerance * max(1e-9, chelators_free[c_idx]):
                    chelators_free[c_idx] = new_c_free
                else:
                    chelators_free[c_idx] = new_c_free
                    converged = False

            # Update metal concentration
            for m_idx in range(num_metals):
                m_tot = metals_total[m_idx]
                if m_tot <= 0:
                    metals_free[m_idx] = 0.0
                    continue

                sum_binding = 0.0
                for c_idx in range(num_chels):
                    x = metal_stoich[m_idx][c_idx]
                    y = ligand_stoich[m_idx][c_idx]
                    if x == 0 and y == 0:
                        continue

                    logk_val = adjusted_logk[m_idx][c_idx]
                    K_val = 10 ** logk_val
                    sum_binding += x * K_val * (metals_free[m_idx] ** (x - 1)) * (chelators_free[c_idx] ** y)

                new_m_free = m_tot / (1.0 + sum_binding)
                diff = abs(new_m_free - metals_free[m_idx])
                if diff < tolerance * max(1e-9, metals_free[m_idx]):
                    metals_free[m_idx] = new_m_free
                else:
                    metals_free[m_idx] = new_m_free
                    converged = False

            if converged:
                break

        return metals_free, chelators_free, iterations
    # end of EquilibriumUtils
    

if __name__ == "__main__":
    # Load protonation and chelator constants from JSON files
    const_utils = ConstantUtils("webmaxc_protonation.json", "webmaxc_chelator.json")

    # Hardcoded inputs, will replace eventually
    chelators = [
        {"name": "EGTA", "concentration_mM": 5.0}
    ]
    metals = [
        {"name": "Ca2", "concentration_mM": 2}
    ]
    ph_value = 7.2
    ionic_strength_value = 0.15
    temp_celsius = 23.0
    
    # Safety check
    if len(chelators) > 1 and len(metals) > 1:
        raise ValueError("This solver cannot have multiple chelators and multiple metals simultaneously.")

    # Create lists to store all chelator-metal pairs
    chelator_names = [chel["name"] for chel in chelators]
    chelator_concs = [chel["concentration_mM"] for chel in chelators]
    metal_names = [metal["name"] for metal in metals]
    metal_concs = [metal["concentration_mM"] for metal in metals]

    # Perform ion adjustments for each chelator-metal pair using a list of dictionaries to store adjusted logKs
    adjusted_logk_matrix = []
    metal_stoich_matrix = []
    ligand_stoich_matrix = []

    # Generate hydrogen concentration and log_f values
    debye_obj = IonAdjuster(
        ph_value, 
        ionic_strength_value, 
        temp_celsius
    )
    
    _, _, _, hydrogen_molar = debye_obj.compute_hydrogen_correction()
    _, _, log_f_std_val = debye_obj.compute_dielectric_and_log_f(
        debye_obj.standard_temperature_celsius, 
        debye_obj.standard_ionic_strength
    )
    _, _, log_f_exp_val = debye_obj.compute_dielectric_and_log_f(
        temp_celsius, 
        ionic_strength_value
    )
    
    temp_kelvin_std = debye_obj.standard_temperature_kelvin
    temp_kelvin_exp = temp_celsius + 273.15    

    for chel_idx, chelator in enumerate(chelators):
        chelator_name = chelator["name"]
        chelator_conc = chelator["concentration_mM"]

        for metal_idx, metal in enumerate(metals):
            metal_name = metal["name"]
            metal_conc = metal["concentration_mM"]

            try:
                constants = const_utils.get_metal_constants(chelator_name, metal_name)
            except ValueError as e:
                print(e)
                continue  

            standard_metal_ligand_log_k = constants["metal_ligand_log_k"]
            standard_metal_ligand_enthalpy = constants["metal_ligand_enthalpy"]
            standard_metal_h_ligand_log_k  = constants["metal_h_ligand_log_k"]
            standard_metal_h_ligand_enthalpy = constants["metal_h_ligand_enthalpy"]
            metal_stoichiometry = constants["metal_stoichiometry"]
            ligand_stoichiometry = constants["ligand_stoichiometry"]

            g_factor_metal = EquilibriumUtils.calc_charge_factor_for_metal_binding(4, 2)  # Example valences
            corrected_metal_ligand_log_k = debye_obj.adjust_log_k_for_ionic_strength_and_temperature(
                standard_metal_ligand_log_k,
                g_factor_metal,
                standard_metal_ligand_enthalpy,
                log_f_std_val,
                log_f_exp_val,
                temp_kelvin_std,
                temp_kelvin_exp
            )

            if standard_metal_h_ligand_log_k != 0:
                corrected_metal_h_ligand_log_k = debye_obj.adjust_log_k_for_ionic_strength_and_temperature(
                    standard_metal_h_ligand_log_k,
                    g_factor_metal,
                    standard_metal_h_ligand_enthalpy,
                    log_f_std_val,
                    log_f_exp_val,
                    temp_kelvin_std,
                    temp_kelvin_exp
                )

            # Generate partial protonation fraction given chelator and pH
            protonation_steps = const_utils.get_protonation_constants(chelator_name)
            total_proton_sites = len(protonation_steps)
            adjusted_proton_log_k_list = []
            for idx, (proto_logk, proto_dh) in enumerate(protonation_steps):
                sites_remaining = total_proton_sites - idx
                g_factor_proton = EquilibriumUtils.calc_charge_factor_for_protonation(sites_remaining)
                corrected_proton_log_k = debye_obj.adjust_log_k_for_ionic_strength_and_temperature(
                    proto_logk,
                    g_factor_proton,
                    proto_dh,
                    log_f_std_val,
                    log_f_exp_val,
                    temp_kelvin_std,
                    temp_kelvin_exp
                )
                adjusted_proton_log_k_list.append(corrected_proton_log_k)
                
            # pH adjustment for binding, combine ML and MLH binding values when available
            apparent_binding_constant = EquilibriumUtils.compute_apparent_binding_constant(
                corrected_metal_ligand_log_k,
                adjusted_proton_log_k_list,
                hydrogen_molar
            )
            
            if standard_metal_h_ligand_log_k != 0:
                apparent_mlh_binding_constant = EquilibriumUtils.compute_apparent_binding_constant(
                    corrected_metal_h_ligand_log_k,
                    adjusted_proton_log_k_list[1:],
                    hydrogen_molar
                )
                
                net_log10_k = math.log10(apparent_binding_constant + apparent_mlh_binding_constant)
            else:
                net_log10_k = math.log10(apparent_binding_constant)

            # Append to matrix
            adjusted_logk_matrix.append(net_log10_k)
            metal_stoich_matrix.append(metal_stoichiometry)
            ligand_stoich_matrix.append(ligand_stoichiometry)

    # Convert stored values to molar
    metals_total = metals_concs = [metal["concentration_mM"] / 1000.0 for metal in metals]  
    chelators_total = [chel["concentration_mM"] / 1000.0 for chel in chelators]  

    # Reshape adjusted_logk_matrix into a list of lists (metals x chelators)
    num_metals = len(metals)
    num_chelators = len(chelators)
    adjusted_logk = []
    metal_stoich = []
    ligand_stoich = []
    for m in range(num_metals):
        logk_row = []
        metal_stoich_row = []
        ligand_stoich_row = []
        for c in range(num_chelators):
            idx = m * num_chelators + c  # Index in the flat list
            logk_row.append(adjusted_logk_matrix[idx])
            metal_stoich_row.append(metal_stoich_matrix[idx])
            ligand_stoich_row.append(ligand_stoich_matrix[idx])
        adjusted_logk.append(logk_row)
        metal_stoich.append(metal_stoich_row)
        ligand_stoich.append(ligand_stoich_row)

    # 13) Solve free concentrations
    metals_free, chelators_free, iterations = EquilibriumUtils.solve_free_metal(
        metals_total,
        chelators_total,
        adjusted_logk,
        metal_stoich,
        ligand_stoich
    )

    # 14) Calculate Kd values for each chelator-metal pair
    kd_values = []
    for m_idx, metal in enumerate(metals):
        for c_idx, chelator in enumerate(chelators):
            net_logk = adjusted_logk[m_idx][c_idx]
            apparent_binding_constant = 10 ** net_logk
            kd = 1.0 / apparent_binding_constant
            kd_values.append({
                "metal": metal["name"],
                "chelator": chelator["name"],
                "Kd": kd
            })

    # 15) Print results
    for m_idx, metal in enumerate(metals):
        for c_idx, chelator in enumerate(chelators):
            print(f"Chelator-Metal Pair: {chelator['name']} - {metal['name']}")
            print(f"  Free Metal [M]       = {metals_free[m_idx]:.6e} M")
            print(f"  Free Chelator [L]    = {chelators_free[c_idx]:.6e} M")
            print(f"  Kd                   = {kd_values[m_idx * num_chelators + c_idx]['Kd']:.6e} M\n")

    print(f"Took {iterations} iterations to converge.")