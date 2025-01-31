import json
import math
import numpy as np
import pandas as pd
from calculations import Calculator, ConstantUtils


class Preprocessing:
    def __init__(
        self,
        proton_json="webmaxc_protonation.json",
        chelator_json="webmaxc_chelator.json"
    ):
        self.proton_json = proton_json
        self.chelator_json = chelator_json
        self.const_utils = ConstantUtils(self.proton_json, self.chelator_json)
        
    @staticmethod
    def get_available_metals_and_chelators(chelator_json):
	    with open(chelator_json, "r") as f:
	        chel_data = json.load(f)

	    chelator_list = []
	    metal_set = set()

	    for block in chel_data["CHELATOR_CONSTANTS"]:
	        c_name = block["chelator"]
	        chelator_list.append(c_name)

	        for metal_name in block["metals"].keys():
	            metal_set.add(metal_name)

	    return sorted(list(metal_set)), sorted(list(set(chelator_list)))

    def compute_equilibrium(self, chelators, metals, pH, ionic_strength, temperature_celsius):
        const_utils = ConstantUtils(self.proton_json, self.chelator_json)
        results = Calculator.mainblock(
            const_utils=const_utils,
            chelators=chelators,
            metals=metals,
            ph_value=pH,
            ionic_strength_value=ionic_strength,
            temp_celsius=temperature_celsius,
        )

        metal_dict = {}
        for m_idx, metal in enumerate(metals):
            metal_name = metal["name"]
            free_conc = results["metals_free"][m_idx]
            if free_conc > 0:
                pM_value = -math.log10(free_conc)
            else:
                pM_value = float("nan")
                free_conc = float("nan")
            metal_dict[metal_name] = (pM_value, free_conc)
        return metal_dict

    def sweep_equilibrium(
        self,
        chelators,
        metals,
        variable="pH",
        start=1.0,
        end=14.0,
        steps=100,
        fixed_pH=7.0,
        fixed_ionic_strength=0.15,
        fixed_temperature=25.0,
        metal_index=0
    ):
        local_metals = [dict(m) for m in metals]
        records = []
        values = np.linspace(start, end, steps)

        for val in values:
            pH = fixed_pH
            ionic_strength = fixed_ionic_strength
            temperature_c = fixed_temperature

            if variable == "pH":
                pH = val
            elif variable == "ionic_strength":
                ionic_strength = val
            elif variable == "temperature":
                temperature_c = val
            elif variable == "metal_conc":
                local_metals[metal_index]["concentration"] = val

            eq_dict = self.compute_equilibrium(
                chelators=chelators,
                metals=local_metals,
                pH=pH,
                ionic_strength=ionic_strength,
                temperature_celsius=temperature_c
            )
            for metal_name, (pM_val, free_conc) in eq_dict.items():
                records.append({
                    "variable_value": val,
                    "Metal": metal_name,
                    "pM": pM_val,
                    "FreeConcentration": free_conc
                })

        df = pd.DataFrame(records)
        df.rename(columns={"variable_value": variable}, inplace=True)
        return df