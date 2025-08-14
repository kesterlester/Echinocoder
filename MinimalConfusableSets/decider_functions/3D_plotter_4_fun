import plotly.graph_objects as go
import plotly.io as pio
from gaussian import generate_gaussian
from sobol_1 import generate_sobol
from electrostatic import generate_electrostatic

vectors_to_plot = generate_electrostatic(40, 3) # need k=3

fig = go.Figure()

for vec in vectors_to_plot:
    fig.add_trace(go.Scatter3d(
        x=[0, vec[0]],
        y=[0, vec[1]],
        z=[0, vec[2]],
        mode='lines+markers',
        marker=dict(size=3),
        line=dict(width=4),
        name=str(vec)
    ))


fig.update_layout(
    scene=dict(
        xaxis=dict(title='X'),
        yaxis=dict(title='Y'),
        zaxis=dict(title='Z')
    ),
    width=800,
    height=800,
    title="3D Vectors Plot", # REMEMBER TO LABEL PLOT
    showlegend=False
)

# Will open web tab
pio.renderers.default = "browser"  
fig.show()
