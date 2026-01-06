import os
import matplotlib.pyplot as plt

def ensure_figures_dir():
    fig_dir = os.path.join(os.path.dirname(__file__), 'figures')
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    return fig_dir

def plot_ssa_trajectories(results):
    fig_dir = ensure_figures_dir()
    scenarios = list(results.keys())

    fig, axs = plt.subplots(1, 3, figsize=(18, 5), sharey=False)
    colors = {'A': 'tab:blue', 'P': 'tab:orange', 'M': 'tab:green', 'Q': 'tab:red'}

    for i, scenario in enumerate(scenarios):
        data = results[scenario]
        t = data['t']

        A = data['A']
        P = data['P']
        M = data['M']
        Q = data['Q']

        for k in range(A.shape[0]):
            axs[i].plot(t, A[k], color=colors['A'], alpha=0.15)
            axs[i].plot(t, P[k], color=colors['P'], alpha=0.15)
            axs[i].plot(t, M[k], color=colors['M'], alpha=0.15)
            axs[i].plot(t, Q[k], color=colors['Q'], alpha=0.15)

        axs[i].plot(t, A.mean(axis=0), color=colors['A'], linewidth=3, label='A (mean)')
        axs[i].plot(t, P.mean(axis=0), color=colors['P'], linewidth=3, label='P (mean)')
        axs[i].plot(t, M.mean(axis=0), color=colors['M'], linewidth=3, label='M (mean)')
        axs[i].plot(t, Q.mean(axis=0), color=colors['Q'], linewidth=3, label='Q (mean)')

        E_h, A_ext = scenario
        axs[i].set_title(f'Time Series at Eh={E_h} and Aext={A_ext}')
        axs[i].set_xlabel('Time')
        if i == 0:
            axs[i].set_ylabel('State value')
        axs[i].legend()

    plt.tight_layout()
    plt.show()
    fig_path = os.path.join(fig_dir, 'ssa_trajectories.png')
    fig.savefig(fig_path)
    print(f"SSA trajectories saved to {fig_path}")