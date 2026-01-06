import argparse
import time
import json
import os

from ode_model import model as ode_model, charts as ode_charts
from ssa_model import model as ssa_model, charts as ssa_charts

def run_odes(base_params, TMAX):
    start_time = time.time()
    ts_results = ode_model.run_time_series(base_params, TMAX)
    es_results = ode_model.run_efflux_sweep(base_params)
    lambda_grid = ode_model.run_invasion(base_params)
    end_time = time.time() - start_time 
    ode_model.analyse_ode_results(base_params)
    ode_charts.plot_time_series(ts_results) 
    ode_charts.plot_phase_bifurcation(base_params, es_results, lambda_grid)
    print(f"Runtime {end_time:.2f} seconds.\n")

def run_ssa(base_params, TMAX, runs):
    start_time = time.time()
    trajectories = ssa_model.run_pdmp_scenarios(base_params, TMAX, runs, ode_model) 
    end_time = time.time() - start_time
    ssa_model.analyse_ssa_results(trajectories) 
    ssa_charts.plot_ssa_trajectories(trajectories)
    print(f"Runtime: {end_time:.2f} seconds.\n")

def main():
    parser = argparse.ArgumentParser(description="Plasmid AMR dynamics model")
    parser.add_argument("--params", type=str, default="parameters/params.json", help="Path to parameter JSON file")
    parser.add_argument("--TMAX", type=int, default=50, help="Maximum time to run models")
    parser.add_argument("--ode", action="store_true", help="Run deterministic ODE simulation")
    parser.add_argument("--ssa", action="store_true", help="Run stochastic SSA simulation")
    parser.add_argument("--runs", type=int, default=50, help="Number of runs for SSA")
    args = parser.parse_args()

    if not os.path.exists(args.params):
        raise FileNotFoundError(f"Parameter file not found: {args.params}")

    with open(args.params, "r") as f:
        base_params = json.load(f)

    if args.ssa:
        run_ssa(base_params, args.TMAX, args.runs)
    else:
        run_odes(base_params, args.TMAX)

if __name__ == "__main__":
    main()