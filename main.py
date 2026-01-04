import json
import argparse
import os
import time

from model.deterministic.invasion import critical_Eh, critical_Aext
from model.deterministic.equilibria import plasmid_free_equilibrium, invasion_eigenvalue
from model.deterministic.odes import run_odes
from model.stochastic.pdmp import run_pdmp_scenarios
from model.stochastic.statistics import extinction_probability, rescue_probability, mean_time_to_resistance
from model.utils.charts import plot_time_series, plot_ssa_trajectories

def get_scenarios():
    scenarios = [
        {
            "name": "No drug + plasmid",
            "A_ext": 0.0,
            "E_h": 0.0
        },
        {
            "name": "Drug + plasmid, no intrinsic efflux",
            "A_ext": 15.0,
            "E_h": 0.0
        },
        {
            "name": "Drug + plasmid + intrinsic efflux",
            "A_ext": 15.0,
            "E_h": 12.0
        }
    ]
    return scenarios

def ode_outcomes(sc, A_star, lambda_P, Eh_crit, Aext_crit):
    print(f"\nScenario: {sc['name']}")
    print(f"Plasmid-free equilibrium (A*, P*, M*, Q*): ({A_star:.1f}, 0, 0, 0)")
    print(f"Plasmid invasion eigenvalue (λP): {lambda_P:.4f}")
    print(f"Critical intrinsic efflux (Eh) = {Eh_crit:.4f}")
    print(f"Critical extracellular drug (Aext) = {Aext_crit:.4f}")

def ode_output(base_params):
    start = time.time()
    ode_results = {}
    scenarios = get_scenarios()
    for sc in scenarios:
        p = base_params.copy()
        p["A_ext"] = sc["A_ext"]
        p["E_h"] = sc["E_h"]

        A_star = plasmid_free_equilibrium(p)
        lambda_P = invasion_eigenvalue(p)
        Eh_crit = critical_Eh(p)
        Aext_crit = critical_Aext(p)

        ode_outcomes(sc, A_star, lambda_P, Eh_crit, Aext_crit)
        ode_results.update(run_odes(p, args.TMAX, sc))

    plot_time_series(ode_results, "ode")
    print(f'\nRuntime: {(time.time() - start):.2f} seconds')

def ssa_outcomes(sc, stats):
    print(f"\nScenario: {sc['name']}")
    print(f"Probability of plasmid extinction (Pext): {stats['P_ext']:.3f}")
    print(f"Probability of resistance rescue (Prescue): {stats['P_rescue']:.3f}")

def ssa_output(base_params):
    start = time.time()
    ssa_results = {}
    scenarios = get_scenarios()

    for sc in scenarios:
        p = base_params.copy()
        p["A_ext"] = sc["A_ext"]
        p["E_h"] = sc["E_h"]

        results = run_pdmp_scenarios(p, args.TMAX, sc, args.steps)
        ssa_results.update(results)

        trajectories = results[sc["name"]]["trajectories"]

        stats = {
            "P_ext": extinction_probability(trajectories),
            "P_rescue": rescue_probability(trajectories),
        }
        ssa_outcomes(sc, stats)

    plot_ssa_trajectories(ssa_results, "ssa")
    print(f'\nRuntime: {(time.time() - start):.2f} seconds')

parser = argparse.ArgumentParser(description="Plasmid efflux dynamics model")
parser.add_argument("--params", type=str, default="parameters/params.json", help="Path to parameter JSON file")
parser.add_argument("--TMAX", type=int, default=50, help="Maximum time to run")
parser.add_argument("--steps", type=int, default=50, help="Number of steps")
args = parser.parse_args()

if not os.path.exists(args.params):
    raise FileNotFoundError(f"Parameter file not found: {args.params}")

with open(args.params, "r") as f:
    base_params = json.load(f)

ode_output(base_params)
ssa_output(base_params)