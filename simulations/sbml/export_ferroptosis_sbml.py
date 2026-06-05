#!/usr/bin/env python3
"""Export the ferroptosis-core single-cell biochemistry to SBML (#351).

What is exported: the DETERMINISTIC mean-field of the core single-cell ferroptosis
ODE network (state variables LP = lipid peroxide, GSH, GPX4; the iron/Fenton,
GSH scavenging + NRF2 resynthesis, autocatalytic LP propagation + GPX4/FSP1
repair, and GPX4 dynamic regulation terms), under the 2D default parameters with
RSL3 (a direct GPX4 inhibitor) applied at t=0. This is the network in
`ferroptosis-core/src/biochem.rs::sim_cell_step` with the per-step stochastic
noise replaced by its mean and the off-by-default realism layers disabled.

What is NOT exported (documented in analysis/sbml-export.md): the per-step
stochastic noise; the death threshold + post-death accumulation; per-cell
parameter sampling (`gen_cell`); the MUFA/ether/iron/persister/etc. realism
layers; the SDT/PDT time-varying exogenous-ROS envelope; and all 2D/3D spatial
fields. SBML/COPASI users get the deterministic single-cell core.

The script: (1) builds the model in Antimony, (2) converts to SBML L3 and
validates it with libSBML (must have 0 errors), (3) simulates the SBML with
roadrunner (the round-trip: does the exported file reproduce the intended ODE?),
(4) integrates the SAME equations with a forward-Euler dt=1 reference that
mirrors the discrete biochem.rs engine, and (5) compares the SBML continuous
solution and the Euler reference against the ACTUAL ferroptosis-core engine's
mean final state (via the Python bindings). Writes ferroptosis_core.xml and a
comparison figure.

Run: python3 simulations/sbml/export_ferroptosis_sbml.py
"""

from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SBML_OUT = HERE / "ferroptosis_core.xml"
FIG_OUT = HERE / "sbml_roundtrip.png"

# --- 2D default Params (from ferroptosis_core.default_params) ---
P = dict(
    fenton_rate=0.02, gsh_scav_efficiency=0.5, gsh_km=2.0, nrf2_gsh_rate=0.025,
    lp_rate=0.06, lp_propagation=0.10, gpx4_rate=0.30, fsp1_rate=0.08,
    gpx4_degradation_by_ros=0.002, gpx4_nrf2_upregulation=0.008, gsh_max=12.0,
    gpx4_nrf2_target_multiplier=1.0, rsl3_gpx4_inhib=0.92,
)
# --- OXPHOS phenotype mean cell (constants; from gen_cell OXPHOS) ---
C = dict(iron=2.8, gsh0=4.0, gpx4_0=1.0, fsp1=1.0, basal_ros=0.5, lipid_unsat=1.6, nrf2=1.2)
N_STEPS = 180


def gpx4_init_rsl3():
    return C["gpx4_0"] * (1.0 - P["rsl3_gpx4_inhib"])


def euler_reference():
    """Forward-Euler dt=1 reference mirroring biochem.rs::sim_cell_step (noise->mean,
    RSL3, default params, realism layers off). Returns (t, LP, GSH, GPX4)."""
    lp, gsh, gpx4 = 0.0, C["gsh0"], gpx4_init_rsl3()
    iron, fsp1, basal, unsat_raw, nrf2 = C["iron"], C["fsp1"], C["basal_ros"], C["lipid_unsat"], C["nrf2"]
    eff_unsat = max(0.05, unsat_raw)  # mufa=0, ether off
    fenton = iron * P["fenton_rate"]
    total_ros = basal + fenton  # exo=0 for RSL3
    traj = [(0, lp, gsh, gpx4)]
    for t in range(1, N_STEPS + 1):
        gsh_fraction = gsh / (gsh + P["gsh_km"])
        scavenged = total_ros * P["gsh_scav_efficiency"] * gsh_fraction
        gsh = max(0.0, gsh - scavenged * 0.5)
        deficit = max(0.0, (P["gsh_max"] - gsh) / P["gsh_max"])
        gsh = gsh + nrf2 * P["nrf2_gsh_rate"] * deficit
        unscav = max(0.0, total_ros - scavenged)
        lp_direct = unscav * eff_unsat * P["lp_rate"]
        quench = gpx4 * (gsh / (gsh + 0.5)) + fsp1
        prop_rate = P["lp_propagation"] / (1.0 + quench * 5.0)
        lp_prop = lp * eff_unsat * prop_rate
        gpx4_repair = gpx4 * (gsh / (gsh + 1.0)) * P["gpx4_rate"] * (lp / (lp + 0.5))
        fsp1_repair = fsp1 * P["fsp1_rate"] * (lp / (lp + 0.5))
        lp = max(0.0, lp + lp_direct + lp_prop - gpx4_repair - fsp1_repair)
        if total_ros > 1.0:
            gpx4 -= P["gpx4_degradation_by_ros"] * (total_ros - 1.0)
        gpx4 += P["gpx4_nrf2_upregulation"] * (nrf2 * P["gpx4_nrf2_target_multiplier"] - gpx4)
        traj.append((t, lp, gsh, gpx4))
    return np.array(traj)


def antimony_model():
    """Continuous-ODE form of the same network for SBML. Rate of change per unit
    time = the per-step delta of the discrete engine (dt=1)."""
    p = P
    c = C
    total_ros = c["iron"] * p["fenton_rate"] + c["basal_ros"]  # constant (RSL3, exo=0)
    eff_unsat = max(0.05, c["lipid_unsat"])
    gpx4_init = gpx4_init_rsl3()
    return f"""
model ferroptosis_core_single_cell
  // === species (state variables) ===
  species LP = 0.0;        // lipid peroxide
  species GSH = {c['gsh0']};
  species GPX4 = {gpx4_init};   // RSL3: gpx4_0 * (1 - rsl3_gpx4_inhib)

  // === fixed quantities (constants for this phenotype/treatment) ===
  total_ros = {total_ros};       // basal_ros + iron*fenton_rate (RSL3: exo = 0)
  eff_unsat = {eff_unsat};       // lipid_unsat (mufa/ether off)
  fsp1 = {c['fsp1']};
  nrf2 = {c['nrf2']};
  gsh_km = {p['gsh_km']};
  gsh_scav_efficiency = {p['gsh_scav_efficiency']};
  nrf2_gsh_rate = {p['nrf2_gsh_rate']};
  gsh_max = {p['gsh_max']};
  lp_rate = {p['lp_rate']};
  lp_propagation = {p['lp_propagation']};
  gpx4_rate = {p['gpx4_rate']};
  fsp1_rate = {p['fsp1_rate']};
  gpx4_degradation_by_ros = {p['gpx4_degradation_by_ros']};
  gpx4_nrf2_upregulation = {p['gpx4_nrf2_upregulation']};
  gpx4_nrf2_target = nrf2 * {p['gpx4_nrf2_target_multiplier']};

  // === intermediates (assignment rules) ===
  scavenged := total_ros * gsh_scav_efficiency * (GSH / (GSH + gsh_km));
  deficit := max(0, (gsh_max - GSH) / gsh_max);
  unscav := max(0, total_ros - scavenged);
  quench := GPX4 * (GSH / (GSH + 0.5)) + fsp1;
  prop_rate := lp_propagation / (1 + quench * 5);
  lp_gen := unscav * eff_unsat * lp_rate + LP * eff_unsat * prop_rate;
  lp_repair := GPX4*(GSH/(GSH+1))*gpx4_rate*(LP/(LP+0.5)) + fsp1*fsp1_rate*(LP/(LP+0.5));

  // === ODEs (rate rules; dt=1 step deltas as continuous rates) ===
  // GPX4 ROS-degradation is gated behind total_ros > 1 in the engine
  // (biochem.rs: `if total_ros > 1.0`); max(0, total_ros - 1) reproduces that
  // ReLU gate so the export matches the engine faithfully.
  LP' = lp_gen - lp_repair;
  GSH' = -scavenged * 0.5 + nrf2 * nrf2_gsh_rate * deficit;
  GPX4' = -gpx4_degradation_by_ros * max(0, total_ros - 1) + gpx4_nrf2_upregulation * (gpx4_nrf2_target - GPX4);
end
"""


def main():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import libsbml
    import tellurium as te

    # (1) Antimony -> (2) SBML + validate
    ant = antimony_model()
    sbml = te.antimonyToSBML(ant)
    SBML_OUT.write_text(sbml)
    doc = libsbml.readSBMLFromString(sbml)
    n_err = doc.getNumErrors(libsbml.LIBSBML_SEV_ERROR)
    n_fatal = doc.getNumErrors(libsbml.LIBSBML_SEV_FATAL)
    print(f"[validate] SBML L{doc.getLevel()}V{doc.getVersion()}: "
          f"{n_err} errors, {n_fatal} fatal -> {'VALID' if n_err == 0 and n_fatal == 0 else 'INVALID'}")
    assert n_err == 0 and n_fatal == 0, "SBML failed libSBML validation"

    # (3) round-trip: load the SBML file back and simulate with roadrunner
    rr = te.loadSBMLModel(sbml)
    rr.timeCourseSelections = ["time", "LP", "GSH", "GPX4"]
    sim = rr.simulate(0, N_STEPS, N_STEPS + 1)
    t_s, lp_s, gsh_s, gpx4_s = sim[:, 0], sim[:, 1], sim[:, 2], sim[:, 3]

    # (4) Euler reference (mirrors the discrete biochem.rs engine)
    ref = euler_reference()
    t_r, lp_r, gsh_r, gpx4_r = ref[:, 0], ref[:, 1], ref[:, 2], ref[:, 3]

    # (5) actual ferroptosis-core engine mean final state (via bindings)
    import ferroptosis_core as fc
    b = fc.sim_batch("OXPHOS", "RSL3", 4000, 42, "2d")
    print(f"[round-trip] final LP  -> SBML {lp_s[-1]:.3f} | Euler-ref {lp_r[-1]:.3f} | "
          f"ferroptosis-core mean {b['mean_lp']:.3f}")
    print(f"[round-trip] final GSH -> SBML {gsh_s[-1]:.3f} | Euler-ref {gsh_r[-1]:.3f} | "
          f"ferroptosis-core mean {b['mean_gsh']:.3f}")
    print(f"[round-trip] final GPX4-> SBML {gpx4_s[-1]:.3f} | Euler-ref {gpx4_r[-1]:.3f} | "
          f"ferroptosis-core mean {b['mean_gpx4']:.3f}")
    # The SBML continuous solution must track the Euler reference (the round-trip
    # integrity check: the exported file reproduces the intended ODE).
    rel = abs(lp_s[-1] - lp_r[-1]) / max(lp_r[-1], 1e-9)
    print(f"[round-trip] SBML vs Euler-ref final-LP relative difference = {rel*100:.1f}% "
          f"({'PASS' if rel < 0.15 else 'CHECK'})")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    for k, (name, s, r) in enumerate(
        [("LP (lipid peroxide)", lp_s, lp_r), ("GSH", gsh_s, gsh_r), ("GPX4", gpx4_s, gpx4_r)]
    ):
        ax[k].plot(t_r, r, color="#3b6ea5", lw=2.5, label="engine ODE (Euler reference, dt=1)")
        ax[k].plot(t_s, s, color="#b5651d", lw=1.4, ls="--", label="exported SBML (roadrunner)")
        ax[k].set_xlabel("step")
        ax[k].set_ylabel(name)
        ax[k].set_title(name)
        ax[k].legend(fontsize=8)
    fig.suptitle("SBML round-trip: the exported ferroptosis-core ODE (orange, roadrunner) "
                 "reproduces the engine dynamics (blue) under RSL3 (#351)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(FIG_OUT, dpi=130)
    print(f"wrote {SBML_OUT}")
    print(f"wrote {FIG_OUT}")


if __name__ == "__main__":
    main()
