"""Rank saved QCD points by a weak-field quartic/binding-energy diagnostic.

Use linear profiles at the existing target energies. The ranking is not a
nonlinear profile prediction. Also plot the actually solved 211 pilot against
its exactly energy-matched linear reference; no joint curve is fabricated.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
import tomllib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import simpson
from scipy.special import lpmv, gammaln

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "src"))
from common import plotstyle

CAMPAIGN = ROOT / "runs/adiabatic/multiphase_two_cloud"


def radial_geometry(positions, spin):
    """Vectorized inversion of the repo's Kerr tortoise-coordinate convention."""
    gap = 2*math.sqrt(1-spin*spin)
    outer_horizon = 1+gap/2
    inner_horizon = spin*spin/outer_horizon
    logarithm = np.where(positions < outer_horizon,
        (positions-outer_horizon+2*inner_horizon/gap*math.log(gap))*gap/(2*outer_horizon),
        np.log(np.maximum(positions-outer_horizon, 1)))
    for _ in range(100):
        distance = np.exp(logarithm)
        radius = outer_horizon+distance
        log_outer = np.logaddexp(logarithm, math.log(gap))
        reconstructed = radius+2*log_outer+2*outer_horizon/gap*(logarithm-log_outer)
        residual = reconstructed-positions
        if np.max(np.abs(residual)) < 1e-10:
            return radius, distance*(distance+gap)/(radius**2+spin**2)
        logarithm -= np.clip(residual*(distance+gap)/(radius**2+spin**2), -2, 2)
    raise ValueError("Tortoise inversion did not converge")


def unit_quartic_diagnostic(profile_path, alpha, spin, m):
    """Evaluate the quartic potential energy on the saved unit-amplitude field.

    All saved linear modes share m: phi=2 Re[S(theta,r) exp(i*m*azimuth)].
    Its exact azimuthal mean is <phi^4>=6*|S|^4. The Kerr potential-energy
    measure then gives dE4/dr*=-(pi*alpha^2/2)*(dr/dr*)*integral rho^2|S|^4 dx.
    """
    with profile_path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    degrees = sorted({int(row["ell"]) for row in rows})
    assert all(int(row["m"]) == m for row in rows)
    polar_nodes, polar_weights = leggauss(24)
    polar = np.array([np.exp(0.5*(math.log((2*l+1)/2)+gammaln(l-m+1)-gammaln(l+m+1)))
                      * lpmv(m, l, polar_nodes) for l in degrees])
    assert np.allclose((polar*polar_weights)@polar.T, np.eye(len(degrees)), atol=1e-13)
    quartic_energy, coarse_energy, peak_field = 0.0, 0.0, 0.0
    for side in ("from_horizon", "from_outer"):
        tables = []
        for degree in degrees:
            selected = sorted((row for row in rows if row["branch"] == side and int(row["ell"]) == degree),
                              key=lambda row: float(row["r_star_over_M"]))
            positions = np.array([float(row["r_star_over_M"]) for row in selected])
            tables.append(np.array([complex(float(row["R_real"]), float(row["R_imag"])) for row in selected]))
        coefficients = np.array(tables)
        radius, jacobian = radial_geometry(positions, spin)
        angular_field = coefficients.T@polar
        peak_field = max(peak_field, 2*np.max(np.abs(angular_field)))
        rho_squared = radius[:,None]**2+spin**2*polar_nodes[None,:]**2
        density = -math.pi*alpha**2/2*jacobian*np.sum(
            polar_weights[None,:]*rho_squared*np.abs(angular_field)**4, axis=1)
        quartic_energy += simpson(density, x=positions)
        coarse_indices = np.unique(np.r_[np.arange(0,len(positions),2),len(positions)-1])
        coarse_energy += simpson(density[coarse_indices], x=positions[coarse_indices])
    assert abs(coarse_energy/quartic_energy-1) < 1e-4
    return quartic_energy, peak_field, abs(coarse_energy/quartic_energy-1)


def plot_pilot(output):
    prefix = CAMPAIGN / "results/qcd_profile_pilots_20261030_v1/qcd211_alpha030_spin099_mass005"
    metadata = tomllib.loads(prefix.with_suffix(".toml").read_text())
    assert metadata["passed"]
    assert math.isclose(metadata["measured_energy"], metadata["linear_comparison_energy"], rel_tol=1e-12)
    with prefix.with_suffix(".csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    curves = {}
    for model in ("linear", "nonlinear"):
        selected = sorted((row for row in rows if row["curve"] == model and row["ell"] == row["m"] == "1"),
                          key=lambda row:(float(row["r_star_over_M"]),row["branch"]))
        curves[model] = np.array([[float(row["r_star_over_M"]),float(row["R_abs"])] for row in selected])
    assert np.array_equal(curves["linear"][:,0],curves["nonlinear"][:,0])
    positions, linear = curves["linear"].T
    omiya = curves["nonlinear"][:,1]
    difference_percent = 100*(omiya-linear)/np.max(linear)
    fig, (upper, lower) = plt.subplots(2,1,figsize=(7,6),sharex=True,
        gridspec_kw={"height_ratios":[2,1],"hspace":0.07})
    upper.plot(positions, omiya, color=plotstyle.PALETTE[1], lw=2, label="Omiya (single cloud)")
    upper.plot(positions, linear, color="0.15", ls="--", lw=1.5, label="Linear")
    upper.set_ylabel(r"$|R_{11}|$")
    upper.legend(frameon=False)
    upper.set_title(r"$M=5M_\odot,\ \alpha=0.3,\ a_*=0.99$"+"\n"+
        rf"Equal energy: $E/(F_a^2M)={metadata['measured_energy']:.5f}$")
    lower.plot(positions,difference_percent,color=plotstyle.PALETTE[1],lw=1.8)
    lower.axhline(0,color="0.5",lw=.8)
    lower.set_ylabel("Omiya − linear\n/ linear peak [%]")
    lower.set_xlabel(r"$r_*/M$")
    for axis in (upper,lower):
        axis.set_xlim(-100,200)
        axis.grid(alpha=.18)
        axis.tick_params(direction="in",top=True,right=True)
    upper.set_ylim(bottom=0)
    for suffix in ("png","pdf"):
        fig.savefig(output/f"available_211_pilot.{suffix}",dpi=180,bbox_inches="tight")
    plt.close(fig)
    return dict(maximum_absolute_difference_over_linear_peak_percent=float(np.max(np.abs(difference_percent))),
        measured_energy=metadata["measured_energy"], joint_curve_computed=False,
        source_prefix=str(prefix.relative_to(ROOT)))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_directory",type=Path)
    args=parser.parse_args()
    population_path=CAMPAIGN/"results/qcd_populations_20261029_v1/population_targets.csv"
    with population_path.open(newline="") as stream:
        populations=list(csv.DictReader(stream))
    cache, records, sources={},[],{population_path,Path(__file__)}
    for row in populations:
        alpha,spin,mass=map(float,(row["alpha"],row["spin"],row["mass_solar"]))
        m=int(row["dominant_m"])
        reference=ROOT/row["primary_reference" if m==1 else "secondary_reference"]
        if reference not in cache:
            data=tomllib.loads(reference.read_text())
            assert data["passed"] and data["normalization"]=="first outer coefficient = 1"
            cache[reference]=(data,unit_quartic_diagnostic(reference.with_suffix(".csv"),alpha,spin,m))
            sources.update((reference,reference.with_suffix(".csv")))
        data,(quartic_unit,field_unit,quadrature_error)=cache[reference]
        energy=float(row["dominant_code_energy"])
        amplitude_squared=energy/data["energy"]
        quartic_energy=quartic_unit*amplitude_squared**2
        binding_fraction=1-data["frequency_real"]/alpha
        score=2*abs(quartic_energy)/(energy*binding_fraction)
        records.append(dict(alpha=alpha,spin=spin,mass_solar=mass,dominant_mode=row["dominant_mode"],
            target_code_energy=energy,quartic_energy_on_linear_field=quartic_energy,
            quartic_over_total_energy=abs(quartic_energy)/energy,
            quartic_to_binding_score=score,
            sampled_linear_field_peak=field_unit*math.sqrt(amplitude_squared),
            radial_quadrature_relative_change=quadrature_error))
    records.sort(key=lambda row:row["quartic_to_binding_score"],reverse=True)
    args.output_directory.mkdir(parents=True,exist_ok=False)
    with (args.output_directory/"visibility_ranking.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(records[0]))
        writer.writeheader();writer.writerows(records)
    pilot=plot_pilot(args.output_directory)
    pilot_prefix=ROOT/pilot["source_prefix"]
    sources.update((pilot_prefix.with_suffix(".csv"),pilot_prefix.with_suffix(".toml"),
                    ROOT/"src/common/plotstyle.py"))
    metadata=dict(cases=len(records),pilot=pilot,
        score="2*abs(E_quartic_on_linear_profile)/(E_linear*(1-omega_R/mu))",
        interpretation="Weak-field interaction/binding diagnostic for ranking; not a solved profile difference or error estimate.",
        joint_population="Other cloud remains in the negligible-classical-population limit at the chosen first-saturation snapshot.",
        normalization="Each linear profile is scaled to its existing population target; plotted pilot references share measured nonlinear energy.",
        source_sha256={str(path.resolve().relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in sources})
    (args.output_directory/"metadata.json").write_text(json.dumps(metadata,indent=2)+"\n")
    (args.output_directory/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    print("Actual pilot:",pilot)
    for row in records[:12]:
        print(row)


if __name__=="__main__":
    main()
