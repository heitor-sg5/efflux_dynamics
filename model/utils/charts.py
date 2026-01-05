import os
import numpy as np
import matplotlib.pyplot as plt

def ensure_dir(model_type):
    out_dir = os.path.join("output", model_type)
    os.makedirs(out_dir, exist_ok=True)
    return out_dir

def plot_time_series(results, model_type):
    fig_dir = ensure_dir(model_type)

    scenario_names = list(results.keys())
    fig, axs = plt.subplots(1, len(scenario_names), figsize=(18, 5))

    for i, name in enumerate(scenario_names):
        sol = results[name]["solution"]
        A, P, M, Q = sol.y

        axs[i].plot(sol.t, A, label="A (drug)")
        axs[i].plot(sol.t, P, label="P (plasmid)")
        axs[i].plot(sol.t, M, label="M (plasmid mRNA)")
        axs[i].plot(sol.t, Q, label="Q (efflux protein)")
        axs[i].set_title(name)
        axs[i].set_xlabel("Time")
        axs[i].set_ylabel("State variables")
        axs[i].legend()

    plt.tight_layout()
    fig_path = os.path.join(fig_dir, "time_series.png")
    fig.savefig(fig_path)
    plt.close(fig)

def plot_ssa_trajectories(results, model_type):
    fig_dir = ensure_dir(model_type)

    scenario_names = list(results.keys())
    fig, axs = plt.subplots(1, len(scenario_names), figsize=(18, 5))

    for i, name in enumerate(scenario_names):
        trajectories = results[name]["trajectories"]

        T_max = max(traj["t"][-1] for traj in trajectories)
        t_common = np.linspace(0, T_max, 500)

        A_all, P_all, M_all, Q_all = [], [], [], []

        for traj in trajectories:
            t = traj["t"]

            A_interp = np.interp(t_common, t, traj["A"])
            P_interp = np.interp(t_common, t, traj["P"])
            M_interp = np.interp(t_common, t, traj["M"])
            Q_interp = np.interp(t_common, t, traj["Q"])

            A_all.append(A_interp)
            P_all.append(P_interp)
            M_all.append(M_interp)
            Q_all.append(Q_interp)

            axs[i].plot(t, traj["A"], alpha=0.3, linewidth=0.8)
            axs[i].plot(t, traj["P"], alpha=0.3, linewidth=0.8)
            axs[i].plot(t, traj["M"], alpha=0.3, linewidth=0.8)
            axs[i].plot(t, traj["Q"], alpha=0.3, linewidth=0.8)

        axs[i].plot(t_common, np.mean(A_all, axis=0), linewidth=3, label="A (mean)")
        axs[i].plot(t_common, np.mean(P_all, axis=0), linewidth=3, label="P (mean)")
        axs[i].plot(t_common, np.mean(M_all, axis=0), linewidth=3, label="M (mean)")
        axs[i].plot(t_common, np.mean(Q_all, axis=0), linewidth=3, label="Q (mean)")

        axs[i].set_title(name)
        axs[i].set_xlabel("Time")
        axs[i].set_ylabel("State variables")
        axs[i].legend()

    plt.tight_layout()
    fig_path = os.path.join(fig_dir, "ssa_trajectories.png")
    fig.savefig(fig_path)
    plt.close(fig)

def plot_growth_and_phase_diagram(Eh_vals, Aext_vals, growth_results, lambda_grid, cP, rP):
    fig_dir = ensure_dir("ode")
    fig, axs = plt.subplots(1, 2, figsize=(14, 5))

    for Aext, data in growth_results.items():
        axs[0].plot( data["Eh"], data["G_star"], label=f"A_ext = {Aext}")

    axs[0].axhline( cP / rP, linestyle="--", color="black", label="cP/rP")
    axs[0].set_xlabel("Intrinsic efflux Eh")
    axs[0].set_ylabel("Growth rate G(A*)")
    axs[0].set_title("Growth rescue by intrinsic efflux")
    axs[0].legend()

    Aext_mesh, Eh_mesh = np.meshgrid(Aext_vals, Eh_vals)

    cmap = axs[1].contourf(Aext_mesh, Eh_mesh, lambda_grid, levels=30, cmap="coolwarm")
    axs[1].contour( Aext_mesh, Eh_mesh, lambda_grid, levels=[0], colors="black", linewidths=2)

    cbar = fig.colorbar(cmap, ax=axs[1])
    cbar.set_label("Invasion eigenvalue λP")

    axs[1].set_xlabel("Extracellular drug Aext")
    axs[1].set_ylabel("Intrinsic efflux Eh")
    axs[1].set_title("Plasmid invasion phase diagram")

    plt.tight_layout()
    plt.savefig(f"{fig_dir}/growth_and_invasion_phase.png")
    plt.close()