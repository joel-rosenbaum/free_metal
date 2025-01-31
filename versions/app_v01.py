import math
import numpy as np
import pandas as pd
import plotly.graph_objects as go

from calculations import Calculator, ConstantUtils


def compute_equilibrium(
    chelators,
    metals,
    pH,
    ionic_strength,
    temperature_celsius,
    proton_json="webmaxc_protonation.json",
    chelator_json="webmaxc_chelator.json"
):
    const_utils = ConstantUtils(proton_json, chelator_json)
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
    chelators,
    metals,
    variable="pH",
    start=1.0,
    end=14.0,
    steps=100,
    fixed_pH=7.0,
    fixed_ionic_strength=0.15,
    fixed_temperature=25.0,
    proton_json="webmaxc_protonation.json",
    chelator_json="webmaxc_chelator.json",
    metal_index=0
):
    """
    Sweeps one variable from 'start' to 'end' in 'steps' increments.
    variable can be:
      - 'pH'
      - 'ionic_strength'
      - 'temperature'
      - 'metal_conc' (varies the metal at position `metal_index` in `metals`).
    Returns a DataFrame: [variable_value, Metal, pM, FreeConcentration].
    """
    local_metals = [dict(m) for m in metals]  # copy
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

        eq_dict = compute_equilibrium(
            chelators=chelators,
            metals=local_metals,
            pH=pH,
            ionic_strength=ionic_strength,
            temperature_celsius=temperature_c,
            proton_json=proton_json,
            chelator_json=chelator_json
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

def x_format(variable):
    return ".2e" if variable == "metal_conc" else ".2f"

def create_dynamic_equilibrium_plot(
    chelators,
    metals,
    proton_json="webmaxc_protonation.json",
    chelator_json="webmaxc_chelator.json"
):
    """
    Single Plotly figure with an updatemenu dropdown for different X-axis sweeps:
      - pH
      - temperature
      - ionic_strength
      - [Metal1 total conc]
      - [Metal2 total conc] ...
    """
    fig = go.Figure()

    sweeps = [
        {
            "label": "pH",
            "variable": "pH",
            "start": 1.0,
            "end": 14.0,
            "steps": 140,
            "fixed_pH": 7.0,
            "fixed_ionic_strength": 0.15,
            "fixed_temperature": 25.0,
            "metal_index": 0,
            "reverse_axes": True
        },
        {
            "label": "temperature",
            "variable": "temperature",
            "start": 10.0,
            "end": 40.0,
            "steps": 30,
            "fixed_pH": 7.0,
            "fixed_ionic_strength": 0.15,
            "fixed_temperature": 25.0,
            "metal_index": 0,
            "reverse_axes": False
        },
        {
            "label": "ionic_strength",
            "variable": "ionic_strength",
            "start": 0.0,
            "end": 0.3,
            "steps": 30,
            "fixed_pH": 7.0,
            "fixed_ionic_strength": 0.15,
            "fixed_temperature": 25.0,
            "metal_index": 0,
            "reverse_axes": False
        }
    ]

    # Add "metal_conc" sweeps for each metal
    for i, metal in enumerate(metals):
        sweeps.append({
            "label": f"[{metal['name']} total]",
            "variable": "metal_conc",
            "start": 1e-6,
            "end": 1e-2,
            "steps": 80,
            "fixed_pH": 7.0,
            "fixed_ionic_strength": 0.15,
            "fixed_temperature": 25.0,
            "metal_index": i,
            "reverse_axes": False
        })

    trace_visibility = []
    trace_index = 0

    # Add traces for each sweep
    for sweep_i, sweep_conf in enumerate(sweeps):
        df_sweep = sweep_equilibrium(
            chelators=chelators,
            metals=metals,
            variable=sweep_conf["variable"],
            start=sweep_conf["start"],
            end=sweep_conf["end"],
            steps=sweep_conf["steps"],
            fixed_pH=sweep_conf["fixed_pH"],
            fixed_ionic_strength=sweep_conf["fixed_ionic_strength"],
            fixed_temperature=sweep_conf["fixed_temperature"],
            metal_index=sweep_conf["metal_index"],
            proton_json=proton_json,
            chelator_json=chelator_json
        )

        local_indices = []
        x_fmt = x_format(sweep_conf["variable"])

        for metal_name, sub_df in df_sweep.groupby("Metal"):
            x_vals = sub_df[sweep_conf["variable"]]
            y_vals = sub_df["pM"]
            free_vals = sub_df["FreeConcentration"]

            # Only first sweep (pH) is visible by default
            initial_visible = (sweep_i == 0)

            # Build the hover template
            hovertemplate = (
                f"{sweep_conf['label']}: %{{x:{x_fmt}}}<br>" 
                f"metal: {metal_name}<br>"
                "pM: %{y:.2f}<br>"
                "[free]: %{customdata:.2e} M<extra></extra>"
            )

            fig.add_trace(
                go.Scatter(
                    x=x_vals,
                    y=y_vals,
                    mode="lines",
                    name=metal_name,
                    customdata=free_vals,
                    hovertemplate=hovertemplate,
                    visible=initial_visible
                )
            )
            local_indices.append(trace_index)
            trace_index += 1

        trace_visibility.append(local_indices)

    num_traces = trace_index

    # Build dropdown
    dropdowns = []
    for sweep_i, sweep_conf in enumerate(sweeps):
        visible_mask = [False] * num_traces
        for idx in trace_visibility[sweep_i]:
            visible_mask[idx] = True

        x_fmt = x_format(sweep_conf["variable"])

        if sweep_conf["reverse_axes"]:
            axis_updates = {
                "xaxis": {
                    "autorange": "reversed",
                    "title": {"text": sweep_conf["label"], "font": {"size": 18}},
                    "tickformat": x_fmt,
                    "tickfont": {"size": 14}
                },
                "yaxis": {
                    "autorange": "reversed",
                    "title": {"text": "pM", "font": {"size": 18}},
                    "tickformat": ".1f",
                    "tickfont": {"size": 14}
                }
            }
        else:
            axis_updates = {
                "xaxis": {
                    "autorange": True,
                    "title": {"text": sweep_conf["label"], "font": {"size": 18}},
                    "tickformat": x_fmt,
                    "tickfont": {"size": 14}
                },
                "yaxis": {
                    "autorange": "reversed",
                    "title": {"text": "pM", "font": {"size": 18}},
                    "tickformat": ".1f",
                    "tickfont": {"size": 14}
                }
            }

        dropdowns.append(
            dict(
                label=sweep_conf["label"],
                method="update",
                args=[
                    {"visible": visible_mask},
                    axis_updates
                ]
            )
        )

    # Initial axis format
    init_x_fmt = x_format(sweeps[0]["variable"])

    fig.update_layout(
        title={
            "text": "Equilibrium Plots",
            "x": 0.10,
            "y": 0.88,
            "xanchor": "left",
            "font_size": 24
        },
        xaxis=dict(
            title=dict(text=sweeps[0]["label"], font=dict(size=18)),
            autorange="reversed" if sweeps[0]["reverse_axes"] else True,
            tickformat=init_x_fmt,
            tickfont=dict(size=14)
        ),
        yaxis=dict(
            title=dict(text="pM", font=dict(size=18)),
            autorange="reversed",
            tickformat=".1f",
            tickfont=dict(size=14)
        ),
        legend=dict(font=dict(size=14)),
        hovermode="closest",
        hoverlabel=dict(
            bgcolor="white",
            font_family="Arial"
        ),
        width=800,
        height=500,
        updatemenus=[
            dict(
                buttons=dropdowns,
                direction="down",
                showactive=True,
                x=1.0,
                y=1.15,
                xanchor="right",
                yanchor="top"
            )
        ]
    )

    return fig


if __name__ == "__main__":
    chelators = [
        {"name": "EGTA", "concentration": 0.005},
    ]
    metals = [
        {"name": "Ca2", "concentration": 0.006},
        {"name": "Zn2", "concentration": 0.001}
    ]

    fig = create_dynamic_equilibrium_plot(chelators, metals)
    fig.show()
