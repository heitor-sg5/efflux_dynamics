import os
import numpy as np
import matplotlib.pyplot as plt

def ensure_figures_dir():
    fig_dir = os.path.join(os.path.dirname(__file__), 'figures')
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    return fig_dir

def plot_time_series(results):
    fig_dir = ensure_figures_dir()
    A_E_values = sorted(results.keys())
    fig, axs = plt.subplots(1, len(A_E_values), figsize=(18, 5), sharey=False)
    for i, (a, e) in enumerate(A_E_values):
        sol = results[(a, e)]
        t = sol.t
        A, P, M, Q = sol.y
        axs[i].plot(t, A, label='A')
        axs[i].plot(t, P, label='P')
        axs[i].plot(t, M, label='M')
        axs[i].plot(t, Q, label='Q')
        axs[i].set_title(f'Time Series at Eh={e}, Aext={a}')
        axs[i].set_xlabel('Time')
        if i == 0:
            axs[i].set_ylabel('State value')
        axs[i].legend()
    plt.tight_layout()
    plt.show()
    fig_path = os.path.join(fig_dir, 'time_series.png')
    fig.savefig(fig_path)
    print(f"Time series figure saved to {fig_path}")

def plot_phase_bifurcation(base_params, results, lambda_grid):
    fig_dir = ensure_figures_dir()
    fig, axs = plt.subplots(1, 2, figsize=(14, 5))

    Eh_vals = np.linspace(0.01, 5.0, 100)
    Aext_vals = np.linspace(0.1, 30.0, 80)

    for Aext, data in results.items():
        axs[0].plot( data["E_h"], data["G_star"], label=f"A_ext = {Aext}")
    axs[0].axhline( base_params['c_P'] / base_params['r_P'], linestyle="--", color="black", label="cP/rP")
    axs[0].set_xlabel("Intrinsic efflux, Eh")
    axs[0].set_ylabel("Growth rate, G(A*)")
    axs[0].set_title("Growth Rescue by Intrinsic Efflux")
    axs[0].legend()

    Aext_mesh, Eh_mesh = np.meshgrid(Aext_vals, Eh_vals)
    cmap = axs[1].contourf(Aext_mesh, Eh_mesh, lambda_grid, levels=30, cmap="coolwarm")
    axs[1].contour( Aext_mesh, Eh_mesh, lambda_grid, levels=[0], colors="black", linewidths=2)
    cbar = fig.colorbar(cmap, ax=axs[1])
    cbar.set_label("Invasion eigenvalue, λP")
    axs[1].set_xlabel("Extracellular drug, Aext")
    axs[1].set_ylabel("Intrinsic efflux, Eh")
    axs[1].set_title("Plasmid Invasion Phase Diagram")

    plt.tight_layout()
    plt.show()
    fig_path = os.path.join(fig_dir, 'phase_bifurcation.png')
    fig.savefig(fig_path)
    print(f"Phase and bifurcation figures saved to {fig_path}")