import dash
from dash import Dash, dcc, html, Input, Output, State, MATCH, ALL, no_update
import plotly.graph_objs as go  
import math
import os

from calculations import Calculator
from processing import Preprocessing
from plotting import Plotting


app = Dash(__name__)
server = app.server 


# Initialize classes
preprocessor = Preprocessing()
const_utils = preprocessor.const_utils
plotter = Plotting()

# Core data and templates
metals_list, chelators_list = preprocessor.get_available_metals_and_chelators("webmaxc_chelator.json")

standard_conditions_dict = {
    "PBS, standard conditions":    {"temp": 25.0, "ionic": 0.1515, "pH": 7.4},
    "Ringer's solution, standard conditions": {"temp": 25.0, "ionic": 0.1572, "pH": 7.0},
    "HBSS, divalent free, incubator": {"temp": 37.0, "ionic": 0.1497, "pH": 7.4},
}
standard_conditions_options = [{"label": k, "value": k} for k in standard_conditions_dict]

header_style = {
    "fontFamily": "Arial, sans-serif",
    "fontSize": "18px",
    "fontWeight": "normal",
    "marginTop": "0px",
    "marginBottom": "10px"
}
body_style = {
    "fontFamily": "Arial, sans-serif",
    "fontSize": "14px",
    "fontWeight": "normal",
    "marginBottom": "6px"
}

# Initialize figure display
fig = plotter.create_dynamic_equilibrium_plot([], [], preprocessor)
fig.update_layout(autosize=True, height=500, margin=dict(l=10, r=0, t=60, b=50))

# Layout
app.layout = html.Div(
    style={"display": "flex"},
    children=[
        # LEFT PANEL
        html.Div(
            style={"width": "30%", "padding": "15px", "display": "flex", "flexDirection": "column"},
            children=[
                # Reagents section
                html.Div(
                    style={"backgroundColor": "#f1f1f1", "padding": "15px", "marginBottom": "20px", "borderRadius": "10px"},
                    children=[
                        html.H3("Reagents", style=header_style),
                        dcc.Dropdown(
                            id="metals-dropdown",
                            options=[{"label": m, "value": m} for m in metals_list],
                            multi=True,
                            placeholder="Select metals",
                            style=body_style
                        ),
                        dcc.Dropdown(
                            id="chelators-dropdown",
                            options=[{"label": c, "value": c} for c in chelators_list],
                            multi=True,
                            placeholder="Select chelators",
                            style=body_style
                        ),
                        html.Div(style={"height": "6px"}),
                        html.Div(id="metals-conc-container", style=body_style),
                        html.Div(id="chelators-conc-container", style=body_style),
                    ]
                ),
                # Conditions section
                html.Div(
                    style={"backgroundColor": "#f1f1f1", "padding": "15px", "marginBottom": "20px", "borderRadius": "10px"},
                    children=[
                        html.H3("Conditions", style=header_style),
                        dcc.Dropdown(
                            id="standard_conditions",
                            options=standard_conditions_options,
                            placeholder="Standard conditions",
                            style=body_style
                        ),
                        html.Div(style={"height": "6px"}),

                        html.Div([
                            html.Label("Temperature: ", style=body_style),
                            dcc.Input(id="temp-input", type="number", placeholder="°C", style=body_style)
                        ]),
                        html.Div([
                            html.Label("Ionic strength: ", style=body_style),
                            dcc.Input(id="ionic-input", type="number", placeholder="M", style=body_style)
                        ]),
                        html.Div([
                            html.Label("pH: ", style=body_style),
                            dcc.Input(id="ph-input", type="number", placeholder="7.0", style=body_style)
                        ]),
                        html.Div(style={"height": "6px"}),
                        html.Button("Calculate", id="calculate-button", style={"fontSize": "18px"}, n_clicks=0)
                    ]
                ),
                html.Div(
                    style={"backgroundColor": "#f1f1f1", "padding": "15px", "borderRadius": "10px"},
                    children=[
                        html.H3("About", style=header_style),
                        html.P([
                            "This tool computes metal-chelator equilibria using calculations derived from the ",
                            html.A(
                                "WEBMAXC calculator",
                                href="https://somapp.ucdmc.ucdavis.edu/pharmacology/bers/maxchelator/webmaxc/webmaxcS.htm", 
                                target="_blank", 
                                style={"color": "#007bff"}
                            ),
                            " originally developed by Chris Patton. Enter your conditions and click CALCULATE to generate ",
                            "free metal concentrations.",
                        ], 
                        style={**body_style, "margin": "0"}) 
                    ]
                )
            ]
        ),
        
        # RIGHT PANEL
        html.Div(
            style={"width": "60%", "padding": "15px"},
            children=[
                # Free species and dissociation constant dropdowns
                html.Div(
                    style={"display": "flex", "justifyContent": "space-between", "marginBottom": "20px"},
                    children=[
                        html.Div([
                            html.H4("Free metal", style=header_style),
                            dcc.Dropdown(id="free-metal-dropdown", style=body_style)
                        ], style={"width": "28%", "backgroundColor": "#f1f1f1", "padding": "15px", "borderRadius": "10px"}),

                        html.Div([
                            html.H4("Free chelator", style=header_style),
                            dcc.Dropdown(id="free-chelator-dropdown", style=body_style)
                        ], style={"width": "28%", "backgroundColor": "#f1f1f1", "padding": "15px", "borderRadius": "10px"}),

                        html.Div([
                            html.H4("Dissociation constant", style=header_style),
                            dcc.Dropdown(id="kd-dropdown", style=body_style)
                        ], style={"width": "28%", "backgroundColor": "#f1f1f1", "padding": "15px", "borderRadius": "10px"}),
                    ]
                ),
                # Graph
                dcc.Graph(figure=fig, id="equilibrium-graph", style={"width": "100%"})
            ]
        )
    ]
)


# Callbacks

@app.callback(
    # Fill pH, ionic, and temp from standard conditions
    Output("temp-input", "value"),
    Output("ionic-input", "value"),
    Output("ph-input", "value"),
    Input("standard_conditions", "value")
)
def fill_in_conditions(selected):
    """
    Whenever a standard condition is selected, update
    the temperature, ionic strength, and pH inputs.
    """
    if not selected:
        return no_update, no_update, no_update

    cond = standard_conditions_dict.get(selected)
    if not cond:
        return no_update, no_update, no_update

    return cond["temp"], cond["ionic"], cond["pH"]


@app.callback(
    Output("metals-conc-container", "children"),
    Output("chelators-conc-container", "children"),
    Input("metals-dropdown", "value"),
    Input("chelators-dropdown", "value"),
)
def update_conc_inputs(metals_selected, chelators_selected):
    """
    Dynamically create the concentration input fields
    for each selected metal and chelator.
    """
    def build_metal_input(metal):
        return html.Div(
            [
                html.Span(f"{metal}: "),
                dcc.Input(
                    id={"type": "metal-conc", "metalName": metal},
                    type="number",
                    placeholder="Conc (M)",
                    style=body_style
                )
            ],
        )

    def build_chelator_input(chel):
        return html.Div(
            [
                html.Span(f"{chel}: "),
                dcc.Input(
                    id={"type": "chelator-conc", "chelatorName": chel},
                    type="number",
                    placeholder="Conc (M)",
                    style=body_style
                )
            ],
        )

    metal_rows = [build_metal_input(m) for m in metals_selected] if metals_selected else []
    chel_rows = [build_chelator_input(c) for c in chelators_selected] if chelators_selected else []
    return metal_rows, chel_rows


@app.callback(
    Output("free-metal-dropdown", "options"),
    Output("free-metal-dropdown", "value"),
    Output("free-chelator-dropdown", "options"),
    Output("free-chelator-dropdown", "value"),
    Output("kd-dropdown", "options"),
    Output("kd-dropdown", "value"),
    Input("calculate-button", "n_clicks"),
    State("ph-input", "value"),
    State("ionic-input", "value"),
    State("temp-input", "value"),
    State("metals-dropdown", "value"),
    State("chelators-dropdown", "value"),
    State({"type": "metal-conc", "metalName": ALL}, "value"),
    State({"type": "metal-conc", "metalName": ALL}, "id"),
    State({"type": "chelator-conc", "chelatorName": ALL}, "value"),
    State({"type": "chelator-conc", "chelatorName": ALL}, "id"),
)
def run_calculation(n_clicks,
                    ph_value, ionic_value, temp_value,
                    metals_selected, chelators_selected,
                    metal_conc_values, metal_conc_ids,
                    chelator_conc_values, chelator_conc_ids):
    """Compute equilibrium and populate the dropdowns for free metals, free chelators, and Kd."""
    if not n_clicks or any(v is None for v in [ph_value, ionic_value, temp_value]) or not all([metals_selected, chelators_selected]):
        return ([], None, [], None, [], None)

    # Build metals/chelators data from user inputs
    metals_data = []
    for conc_val, conc_id in zip(metal_conc_values, metal_conc_ids):
        metal_name = conc_id["metalName"]
        if metal_name in metals_selected:
            c = float(conc_val) if conc_val is not None else 0.0
            metals_data.append({"name": metal_name, "concentration": c})

    chelators_data = []
    for conc_val, conc_id in zip(chelator_conc_values, chelator_conc_ids):
        chel_name = conc_id["chelatorName"]
        if chel_name in chelators_selected:
            c = float(conc_val) if conc_val is not None else 0.0
            chelators_data.append({"name": chel_name, "concentration": c})

    if not metals_data or not chelators_data:
        return ([], None, [], None, [], None)

    # Perform the calculation
    result = Calculator.standard_output(
        const_utils=const_utils,
        metals=metals_data,
        chelators=chelators_data,
        ph_value=ph_value,
        ionic_strength_value=ionic_value,
        temp_celsius=temp_value
    )

    # Prepare dropdown options
    species_data = result["species_data"]  # [{"name": ..., "free": ..., "total": ...}, ...]
    complex_data = result["complex_data"]  # [{"complex": ..., "Kd": ..., ...}, ...]

    metal_options = []
    chelator_options = []
    for sp in species_data:
        label_str = f"{sp['name']}: {sp['free']:.2e} M"
        if sp["name"] in metals_selected:
            metal_options.append({"label": label_str, "value": sp["name"]})
        else:
            chelator_options.append({"label": label_str, "value": sp["name"]})

    kd_options = []
    for cd in complex_data:
        kd_lbl = f"{cd['complex']}: {cd['Kd']:.2e} M"
        kd_options.append({"label": kd_lbl, "value": cd["complex"]})

    # Select first by default
    metal_value    = metal_options[0]["value"] if metal_options else None
    chelator_value = chelator_options[0]["value"] if chelator_options else None
    kd_value       = kd_options[0]["value"] if kd_options else None

    return (
        metal_options,   metal_value,
        chelator_options, chelator_value,
        kd_options,      kd_value
    )


@app.callback(
    Output("equilibrium-graph", "figure"),
    Input("calculate-button", "n_clicks"),
    State("ph-input", "value"),
    State("ionic-input", "value"),
    State("temp-input", "value"),
    State("metals-dropdown", "value"),
    State("chelators-dropdown", "value"),
    State({"type": "metal-conc", "metalName": ALL}, "value"),
    State({"type": "metal-conc", "metalName": ALL}, "id"),
    State({"type": "chelator-conc", "chelatorName": ALL}, "value"),
    State({"type": "chelator-conc", "chelatorName": ALL}, "id"),
    prevent_initial_call=True
)
def update_equilibrium_plot(n_clicks,
                            ph_value, ionic_value, temp_value,
                            metals_selected, chelators_selected,
                            metal_conc_values, metal_conc_ids,
                            chelator_conc_values, chelator_conc_ids):
    """After compute, re-generate the graph with updated data."""
    if not n_clicks:
        raise dash.exceptions.PreventUpdate

    # If missing fields or nothing selected, return an empty figure
    if (ph_value is None or ionic_value is None or temp_value is None 
        or not metals_selected or not chelators_selected):
        return go.Figure()

    # Build metals/chelators data
    metals_data = [
        {"name": idx["metalName"], "concentration": float(val) if val else 0.0}
        for val, idx in zip(metal_conc_values, metal_conc_ids)
        if idx["metalName"] in metals_selected
    ]
    chelators_data = [
        {"name": idx["chelatorName"], "concentration": float(val) if val else 0.0}
        for val, idx in zip(chelator_conc_values, chelator_conc_ids)
        if idx["chelatorName"] in chelators_selected
    ]

    # Generate the updated figure
    fig = plotter.create_dynamic_equilibrium_plot(
        chelators=chelators_data,
        metals=metals_data,
        preprocessor=preprocessor,
        user_pH=ph_value,
        user_ionic=ionic_value,
        user_temp=temp_value
    )
    return fig


if __name__ == "__main__":
    port = int(os.getenv("PORT", 8050))  
    app.run_server(debug=False, host="0.0.0.0", port=port)

server = app.server
    

