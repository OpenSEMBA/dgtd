# Holland wire agreement with FDTD

The 12-rank DGTD run of `3D_Holland1981` and the FDTD case `holland1981` agree on the incident wave and on the wire geometry. They do not agree on the current. This note records that comparison and the solver changes that could close the gap. The incident-wave path stays as it is.

Numbers below are the finished runs: DGTD `exports/SimulationData/mpi-12/3D_Holland1981`, FDTD `holland1981.fdtd_mid_point_Wz_11_11_12_s8_tm.dat`.

## What already agrees

The off-wire probe at `(1.5, 1.1, 1.2)` peaks at $E_z = 0.236$. The file [holland.exc](../../testData/maxwellInputs/3D_Holland1981/holland.exc) peaks at $0.238$. $|H|$ peaks at the same value as $|E|$. That is this solver's normalization, $Z_0 = 1$ inside the volume update.

The wire is the FDTD segment. It runs from `(1.1, 1.1, 0.7)` to `(1.1, 1.1, 1.7)`, length $1\,\mathrm{m}$. The probe is the midpoint `(1.1, 1.1, 1.2)`.

The current peaks at $6.66\,\mathrm{ns}$ here and at $6.12\,\mathrm{ns}$ in FDTD. The file itself peaks at $2.97\,\mathrm{ns}$. Travel from the plane $x = 0$ to the wire is $1.1/c = 3.67\,\mathrm{ns}$, and $2.97 + 3.67 = 6.64\,\mathrm{ns}$. FDTD starts that clock on the upstream corner of its total-field box, which is why its peak is about $0.5\,\mathrm{ns}$ earlier. `magnitude.file` keeps time zero on $\mathbf{x}\cdot\hat{\mathbf{d}} = 0$. That shift is not part of the fix list.

FDTD prints a warning that a wire ending is assumed PEC. That print does not change the terminal. Both ends are logged as embedded in air. Both models intend open ends.

## What does not agree

| Quantity | FDTD | DGTD |
|---|---|---|
| Mid-wire peak current | $4.34\times 10^{-4}$ | $5.30\times 10^{-3}$ |
| Ratio | | about $12$ |
| Correlation, $0$ to $30\,\mathrm{ns}$ | | $0.18$ |
| Sign changes after the peak | about every $4.1\,\mathrm{ns}$ (near $120\,\mathrm{MHz}$) | first zero at $15.9\,\mathrm{ns}$, still negative at $30\,\mathrm{ns}$ |

A $1\,\mathrm{m}$ open dipole rings at $c/(2L) = 150\,\mathrm{MHz}$, a period of $6.7\,\mathrm{ns}$. FDTD is near that. This solver is not.

FDTD's tip at $z = 1.7$ lies in the PML. A host tetrahedron in PML is now accepted when `kappa_max` is 1: the stretch factor on $E$ is then 1, the wire still uses physical arc length and physical $E$, and the current is added on the electric rows after the absorber. `kappa_max` above 1 is rejected. The Holland mesh still places that tip in vacuum, so this acceptance does not by itself move the comparison.

## Unit map

FDTD steps the segment current in seconds. In `wiresconstantes` the electric coefficient is $c_{te2} = \Delta t / L$, with Boutayeb's $L$ already carrying $\mu_0$:

$$
\partial_t I = E / L,\qquad L = \mu_0\, L_\mathrm{hat}
$$

$L_\mathrm{hat}$ is the dimensionless bracket over $4\pi$. For equal transverse steps $\Delta$ and radius $a$,

$$
L_\mathrm{hat} = \frac{1}{4\pi}\left(\log\frac{2\Delta^2}{4a^2} + \frac{\pi}{2} + \pi\frac{a^2}{\Delta^2} - 3\right)
$$

and, when $a < 0.3\Delta$, subtract $0.57/(4\pi)$. That bracket was the inductance used for the comparison below. The live wire uses Edelvik's tube: $r_0=1.7\Delta$ and $L_\mathrm{hat}=(1/(2\pi))\log((r_0+a)/(2a))$, in `edelvikNormalizedInductance`.

This solver steps in light-metres, $t' = c\, t$, and divides by $L_\mathrm{hat}$ rather than by $\mu_0 L_\mathrm{hat}$:

$$
\partial_{t'} I = E / L_\mathrm{hat}
$$

The SI rate is then $c/L_\mathrm{hat}$, which is $Z_0 \approx 377$ times FDTD's rate. If the wire were a lumped inductor, our current would be about $377$ times the FDTD current. The Holland peak is only about $12$ times larger ($377/12 \approx 31$). Radiation and the return into $E$ set the current, not $L$ alone.

Dividing the probe by $Z_0$ would label the number as amperes and would make it about $30$ times smaller than FDTD. That division is not the fix.

The volume return in that comparison was $-M_E^{-1}$ times the line integral of $I$ on the host tetrahedron. The live return is the mass-inverted transpose of Edelvik's tube weight. FDTD subtracts $I$ from the single grid edge through

$$
c_{te5} = \Delta t / (\varepsilon\, \Delta_1 \Delta_2)
$$

## Cell size in $L_\mathrm{hat}$

FDTD's two transverse steps on this case are $\Delta = 0.1\,\mathrm{m}$. With $a = 0.02\,\mathrm{m}$,

$$
L_\mathrm{hat}(0.1\,\mathrm{m}) = 0.0519
$$

Host tetrahedra along this wire have mean edges from $0.113\,\mathrm{m}$ to $0.148\,\mathrm{m}$, average $0.127\,\mathrm{m}$:

$$
L_\mathrm{hat}(0.127\,\mathrm{m}) = 0.0866
$$

That is $1.7$ times the FDTD inductance. It moves the current, and it does not explain a factor of $12$.

## Fixes, in order

### 1. Report the scale, and leave the run alone

The probe file and the ParaView field `current` are the light-metre unknown, not amperes. A later column $I_\mathrm{SI} = I / Z_0$ is only a label. It is not expected to match FDTD. Do not change `rhs` for this step.

### 2. Optional transverse length

FDTD passes the two grid steps into the log. The comparison above used the host mean edge in the Boutayeb bracket. For $a = 0.02\,\mathrm{m}$ that bracket moves from $0.052$ at $0.1\,\mathrm{m}$ to $0.087$ at $0.127\,\mathrm{m}$.

The JSON field is `transverseStep`, in metres. Omit it and the host mean edge is $\Delta$ in the live Edelvik tube. `3D_Wire_Box` no longer uses the Boutayeb value of $L_\mathrm{hat}$. On this Holland mesh, `transverseStep` $0.1$ with `time_step` $0.004$ drove the current to overflow at $5.3\,\mathrm{ns}$, so the case leaves the field unset.

### 3. Check the electrical length of the open wire

The open end sets the incoming characteristic so that $\hat{I} = 0$. FDTD's open node uses the image $I_+ = -I_-$. Those are the same intention. They are not the same stencil.

The check is a $1\,\mathrm{m}$ wire in a plain box, no PML and no TFSF, driven by a short pulse. The current spectrum should peak near $150\,\mathrm{MHz}$. If it peaks far below that, the bug is the line mass, the stiffness, or the end flux in `rhs`, not the Holland mesh and not the excitation file.

### 4. Scale the drive and the return together

Only after fix 3 shows a $150\,\mathrm{MHz}$ wire, and only if the Holland amplitude is still high.

The unknown can be redefined as amperes. Then the electric drive and the return into $E$ are both multiplied by $1/Z_0$. One side alone changes the radiation resistance. The lumped estimate says a drive-only factor of $1/Z_0$ would land about $30$ times under FDTD, so it is forbidden as a one-line change.

### 5. Narrow the return, only if 2–4 are not enough

FDTD's source is one grid edge. Ours is the whole host tetrahedron, which spreads $I$ and weakens the near field. A narrower support would raise the local $E$ reaction and lower the current. That edits the wire source added after item 9 of `Mult`, which is a frozen contract. It waits for a separate decision. This note does not implement it.

## Out of scope

`magnitude.file` already places time zero on $\mathbf{x}\cdot\hat{\mathbf{d}} = 0$. Moving that origin by half a nanosecond to imitate FDTD's upstream corner does not address the factor of $12$ or the missing $150\,\mathrm{MHz}$ ring.
