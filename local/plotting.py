import math
import plotly.graph_objects as go
from calculations import Calculator

class Plotting:
    def __init__(self):
        pass

    @staticmethod
    def x_format(variable):
        return ".2e" if variable == "metal_conc" else ".2f"

    def _create_empty_figure(self):
        """Return an empty figure with a minimal, responsive layout."""
        fig = go.Figure()
        fig.update_layout(
            autosize=True,
            margin=dict(l=10, r=10, t=10, b=10),
            xaxis=dict(
                showgrid=False,
                zeroline=False,
                showticklabels=False
            ),
            yaxis=dict(
                showgrid=False,
                zeroline=False,
                showticklabels=False
            ),
            plot_bgcolor="#f5f5f5"
        )
        return fig

    def create_dynamic_equilibrium_plot(
        self,
        chelators,
        metals,
        preprocessor,
        user_pH=7.0,
        user_ionic=0.15,
        user_temp=25.0
    ):
        """
        Create a dynamic equilibrium plot with "sweeps" over pH, temperature, ionic_strength,
        and metal concentration. Uses a responsive (autosize) layout.
        """
        # If no data, just return an empty figure
        if not metals or not chelators:
            return self._create_empty_figure()

        fig = go.Figure()

        # Sweep (range) configurations
        sweeps = [
            {
                "label": "pH",
                "variable": "pH",
                "start": 4.0,
                "end": 10.0,
                "steps": 300,
                "fixed_pH": user_pH,
                "fixed_ionic_strength": user_ionic,
                "fixed_temperature": user_temp,
                "metal_index": 0,
                "reverse_axes": True
            },
            {
                "label": "temperature",
                "variable": "temperature",
                "start": 0,
                "end": 50.0,
                "steps": 100,
                "fixed_pH": user_pH,
                "fixed_ionic_strength": user_ionic,
                "fixed_temperature": user_temp,
                "metal_index": 0,
                "reverse_axes": False
            },
            {
                "label": "ionic_strength",
                "variable": "ionic_strength",
                "start": 0.0,
                "end": 1.0,
                "steps": 100,
                "fixed_pH": user_pH,
                "fixed_ionic_strength": user_ionic,
                "fixed_temperature": user_temp,
                "metal_index": 0,
                "reverse_axes": False
            }
        ]

        # Dynamically add sweeps for the total concentration of each metal
        for i, metal in enumerate(metals):
            sweeps.append({
                "label": f"[{metal['name']} total]",
                "variable": "metal_conc",
                "start": 1e-6,
                "end": 1e-2,
                "steps": 200,
                "fixed_pH": user_pH,
                "fixed_ionic_strength": user_ionic,
                "fixed_temperature": user_temp,
                "metal_index": i,
                "reverse_axes": False
            })

        # Build traces for each sweep   
        trace_visibility = []
        trace_index = 0

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

                # Only the first sweep is visible initially
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

        # Update menu items for each sweep
        dropdowns = []
        for sweep_i, sweep_conf in enumerate(sweeps):
            visible_mask = [False] * num_traces
            for idx in trace_visibility[sweep_i]:
                visible_mask[idx] = True

            x_fmt = self.x_format(sweep_conf["variable"])

            # Handling for reverse (pH) or normal axis
            if sweep_conf["reverse_axes"]:
                axis_updates = {
                    "xaxis": {
                        "autorange": "reversed",
                        "title": {"text": sweep_conf["label"], "font": {"size": 16}},
                        "tickformat": x_fmt,
                        "tickfont": {"size": 12}
                    },
                    "yaxis": {
                        "autorange": "reversed",
                        "title": {"text": "pM", "font": {"size": 16}},
                        "tickformat": ".1f",
                        "tickfont": {"size": 12}
                    }
                }
            else:
                axis_updates = {
                    "xaxis": {
                        "autorange": True,
                        "title": {"text": sweep_conf["label"], "font": {"size": 16}},
                        "tickformat": x_fmt,
                        "tickfont": {"size": 12}
                    },
                    "yaxis": {
                        "autorange": "reversed",
                        "title": {"text": "pM", "font": {"size": 16}},
                        "tickformat": ".1f",
                        "tickfont": {"size": 12}
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

        # Default to the first sweep (pH)
        init_x_fmt = self.x_format(sweeps[0]["variable"])

        # Position the dropdown
        dropdown_pos = 0.35 if len(metals) < 2 else 0.39

        fig.update_layout(
            autosize=True,
            margin=dict(l=10, r=10, t=50, b=50),
            template="seaborn",
            plot_bgcolor="#f5f5f5",

            # Axes
            xaxis=dict(
                title=dict(text=sweeps[0]["label"], font=dict(size=16)),
                autorange="reversed" if sweeps[0]["reverse_axes"] else True,
                tickformat=init_x_fmt,
                tickfont=dict(size=12)
            ),
            yaxis=dict(
                title=dict(text="pM", font=dict(size=16)),
                autorange="reversed",
                tickformat=".1f",
                tickfont=dict(size=12)
            ),

            # Title above the dropdown
            annotations=[
                dict(
                    x=0.01,
                    y=1.12,
                    xref="paper",
                    yref="paper",
                    text="Free metal as a function of",
                    showarrow=False,
                    font=dict(size=18),
                    align="left"
                )
            ],

            # Dropdown
            updatemenus=[
                dict(
                    buttons=dropdowns,
                    direction="down",
                    showactive=True,
                    x=dropdown_pos,
                    y=1.13,
                    xanchor="left",
                    yanchor="top"
                )
            ],

            legend=dict(font=dict(size=12)),
            hovermode="closest",
            hoverlabel=dict(bgcolor="white", font_family="Arial")
        )

        return fig
