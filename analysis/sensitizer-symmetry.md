# The sonosensitizer was exempt from the barriers RSL3 faced (#SENS-SYM)

## The defect

`simulations/sim-tme/src/main.rs` applied the pH ion-trapping penalty inside a
branch gated on the treatment arm:

```rust
if tx == Treatment::RSL3 {
    // ... reduce effective GPX4 inhibition in acidic zones ...
}
```

There was no corresponding branch for SDT or PDT, so the sensitizer sat at
uniform concentration everywhere while the pharmacologic drug was penalised by
the acidity gradient. Manuscript §7.4 justified this:

> SDT is barely affected (+0.8%, from 139,640 to 140,693 ferroptosis kills)
> because ultrasound energy delivery is pH-independent — **there is no drug to
> trap**.

That sentence is a category error. SDT and PDT both require a *sensitizer*,
which is a systemically dosed small molecule; the manuscript's own lead clinical
exemplar, SONALA-001, is 5-ALA. The ultrasound or light *field* is
pH-independent. The molecule that converts that field into ROS is not.

So part of the reported cross-modality resistance gap was produced by which
agent had been given delivery physics, not by biology.

## The fix

`PhConfig` gains `sensitizer_ion_trap_sensitivity`. When it is non-zero, the
consumer scales each cell's `exo_ros_peak` by the same clamped availability
factor RSL3's GPX4 correction uses:

```
availability = clamp(1 - sensitizer_ion_trap_sensitivity * (ph_edge - local_ph), 0.3, 1.0)
```

Setting it equal to `ion_trap_sensitivity` (0.4) gives the sensitizer exactly
the barrier the drug faces. `0.0` is the default and reproduces the historical
behaviour byte-identically, so the production matrix is unchanged.

Reproduce with:

```bash
cargo build --release -p sim-tme
./target/release/sim-tme                                 # asymmetric (historical)
FERRO_SENSITIZER_ION_TRAP=0.4 ./target/release/sim-tme   # symmetric
```

## Result

pH-gradient condition, λ=120 μm, immune on, seed 42:

| | asymmetric (published) | symmetric | change |
|---|---|---|---|
| RSL3 ferroptosis kills | 77 | 77 | unchanged |
| RSL3 overall kill | 0.0% | 0.0% | unchanged |
| **SDT ferroptosis kills** | **140,693** | **90,403** | **−35.7%** |
| **SDT overall kill** | **88.8%** | **57.7%** | **−31.1 pts** |
| SDT normoxic kill | 95.6% | 91.1% | −4.5 pts |
| **SDT hypoxic kill** | **87.7%** | **53.1%** | **−34.6 pts** |

The manuscript's "+0.8%, barely affected" becomes **−36%** once the sensitizer
is subject to the same acidity barrier as the drug.

## What survives, and what does not

**The direction survives.** RSL3 still collapses to ~0% in the acidic gradient
while SDT retains 57.7%. Physical ROS delivery remains less vulnerable to
ion-trapping than pharmacologic ferroptosis induction, because only part of the
SDT chain is a small molecule while essentially all of the RSL3 chain is.

**The magnitude does not.** "Barely affected" was an artifact. A modality that
loses a third of its kill and half its hypoxic-zone kill to tumour acidity is
not unaffected by that barrier, and §7.4's supporting sentence — "there is no
drug to trap" — is wrong as written.

This is the outcome the manuscript's own abstract says it is aiming for:
"the defensible contributions are directional and methodological, not
magnitudes". Here that framing is load-bearing rather than decorative.

## Scope and limits

- Only the **pH** leg was asymmetric in `sim-tme`. The CAF/stromal boost is
  applied to all arms with no treatment gate, and RSL3 receives no spatial
  penetration gradient in this binary either, so the earlier audit's claim that
  Krogh penetration and CAF supply were also RSL3-only does not hold for this
  code path. The pH leg is the real one, and it is now fixed.
- `sensitizer_ion_trap_sensitivity = 0.4` assumes the sensitizer has the same
  ion-trapping behaviour as RSL3. That is a *symmetry assumption*, not a
  measurement: 5-ALA/PpIX ionisation differs from a chloroacetamide's, and
  `ion_trap_sensitivity` is itself flagged in §7.4 as "the most uncertain
  parameter in the entire model". The defensible statement is that the true
  sensitizer penalty is somewhere between 0 and the drug's, and the published
  number assumed the zero end without saying so.
- Single seed (42), as with every number in this suite. The seed-replication
  gap is separate and still open.
- `sim-tme-3d` carries the same `tx == Treatment::RSL3` pattern at
  `main.rs:1363` and has **not** been changed here.

## Consequence for the manuscript

§7.4 and the §7.5 cumulative-asymmetry argument overstate the pH leg. §7.4 has
been corrected to report both arms. The hypoxia leg (§7.1) is a separate,
already-contested question and is untouched by this fix.
