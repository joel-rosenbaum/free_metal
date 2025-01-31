import math
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import dash
from dash import dcc, html

from calculations import Calculator, ConstantUtils


#######################
# 1. Preprocessing Class
#######################
class Preprocessing:
    def __init__(
        self,
        proton_json="webmaxc_protonation.json",
        chelator_json="webmaxc_chelator.json"
    ):
        self.proton_json = proton_json
        self.chelator_json = chelator_json

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


##################
# 2. Plotting Class
##################
class Plotting:
    def __init__(self):
        pass
        
    @staticmethod
    def x_format(variable):
        return ".2e" if variable == "metal_conc" else ".2f"

    def create_dynamic_equilibrium_plot(self, chelators, metals, preprocessor):
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

        # Dynamically add "metal_conc" for each metal
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

        # Build each sweep's data/traces
        for sweep_i, sweep_conf in enumerate(sweeps):
            df_sweep = preprocessor.sweep_equilibrium(
                chelators=chelators,
                metals=metals,
                variable=sweep_conf["variable"],
                start=sweep_conf["start"],
                end=sweep_conf["end"],
                steps=sweep_conf["steps"],
                fixed_pH=sweep_conf["fixed_pH"],
                fixed_ionic_strength=sweep_conf["fixed_ionic_strength"],
                fixed_temperature=sweep_conf["fixed_temperature"],
                metal_index=sweep_conf["metal_index"]
            )

            local_indices = []
            x_fmt = self.x_format(sweep_conf["variable"])

            for metal_name, sub_df in df_sweep.groupby("Metal"):
                x_vals = sub_df[sweep_conf["variable"]]
                y_vals = sub_df["pM"]
                free_vals = sub_df["FreeConcentration"]

                initial_visible = (sweep_i == 0)
                
                display_label = sweep_conf["label"]
                if display_label.endswith("total]"):
                    display_label = display_label[:-7] + "] (total)"
                
                hovertemplate = (
                    f"{display_label}: %{{x:{x_fmt}}}<br>"
                    f"[{metal_name}] (pM): %{{y:.2f}}<br>"
                    f"[{metal_name}] (free): %{{customdata:.2e}}<extra></extra>"
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

            x_fmt = self.x_format(sweep_conf["variable"])

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

        # Initial x-axis format
        init_x_fmt = self.x_format(sweeps[0]["variable"])

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
            template="seaborn",      
            plot_bgcolor="#f5f5f5",
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


###################################
# 3. Build the Dash Application
###################################
app = dash.Dash(__name__)

preprocessor = Preprocessing()
plotter = Plotting()

# Example data
chelators = [
    {"name": "EGTA", "concentration": 0.005},
]
metals = [
    {"name": "Ca2", "concentration": 0.006},
    {"name": "Zn2", "concentration": 0.001}
]

fig = plotter.create_dynamic_equilibrium_plot(chelators, metals, preprocessor)

###############################
# 4. The Two-Column Layout
###############################
app.layout = html.Div(
    style={
        "display": "flex",
        "flexDirection": "row",
        "minHeight": "100vh"
    },
    children=[
        # LEFT SIDEBAR
        html.Div(
            style={
                "width": "30%",
                "padding": "20px",
                "display": "flex",
                "flexDirection": "column",
                "alignItems": "flex-start",
            },
            children=[
                # -- FIRST RECTANGLE (Metals & Chelators) --
                html.Div(
                    style={
                        "backgroundColor": "#f1f1f1",
                        "borderRadius": "10px",
                        "padding": "15px",
                        "width": "100%",        # fill parent's width
                        "marginBottom": "20px"  # spacing from second rectangle
                    },
                    children=[
                        html.H3("Metals and Chelators"),
                        # Multi-select dropdown for metals
                        dcc.Dropdown(
                            id="metals-dropdown",
                            options=[
                                {"label": "Ca2", "value": "Ca2"},
                                {"label": "Zn2", "value": "Zn2"},
                                {"label": "Mg2", "value": "Mg2"},
                                # etc. 
                            ],
                            multi=True,
                            placeholder="Select metals"
                        ),
                        html.Br(),
                        # Multi-select dropdown for chelators
                        dcc.Dropdown(
                            id="chelators-dropdown",
                            options=[
                                {"label": "EGTA", "value": "EGTA"},
                                {"label": "EDTA", "value": "EDTA"},
                                {"label": "BAPTA", "value": "BAPTA"},
                                # etc.
                            ],
                            multi=True,
                            placeholder="Select chelators"
                        ),
                        html.Br(),
                        # Reserve up to 5 lines for user-specified concentrations
                        # (no callback logic yet — just placeholders)
                        html.Div([
                            html.Div([
                                html.Span("Ca2: "),
                                dcc.Input(type="number", placeholder="Conc (M)", style={"width": "90px"})
                            ], style={"marginBottom": "5px"}),
                            html.Div([
                                html.Span("Zn2: "),
                                dcc.Input(type="number", placeholder="Conc (M)", style={"width": "90px"})
                            ], style={"marginBottom": "5px"}),
                            html.Div([
                                html.Span("EGTA: "),
                                dcc.Input(type="number", placeholder="Conc (M)", style={"width": "90px"})
                            ]),
                            # etc. up to 5
                        ])
                    ]
                ),

                # -- SECOND RECTANGLE (Conditions) --
                html.Div(
                    style={
                        "backgroundColor": "#f1f1f1",
                        "borderRadius": "10px",
                        "padding": "15px",
                        "width": "100%"
                    },
                    children=[
                        html.H3("Conditions"),
                        # Temperature
                        html.Div([
                            html.Label("Temperature (°C):", style={"marginRight": "5px"}),
                            dcc.Input(type="number", placeholder="25.0", style={"width": "90px"})
                        ], style={"marginBottom": "10px"}),
                        
                        # Ionic Strength
                        html.Div([
                            html.Label("Ionic Strength:", style={"marginRight": "5px"}),
                            dcc.Input(type="number", placeholder="0.15", style={"width": "90px"})
                        ], style={"marginBottom": "10px"}),

                        # pH
                        html.Div([
                            html.Label("pH:", style={"marginRight": "5px"}),
                            dcc.Input(type="number", placeholder="7.0", style={"width": "90px"})
                        ]),
                    ]
                )
            ]
        ),

        # RIGHT MAIN PANEL
        html.Div(
            style={
                "width": "70%",
                "padding": "20px",
            },
            children=[
                # The dynamic figure
                dcc.Graph(
                    figure=fig,
                    id="equilibrium-graph"
                )
            ]
        )
    ]
)


if __name__ == "__main__":
    app.run_server(debug=True, use_reloader=False, host="127.0.0.1")
