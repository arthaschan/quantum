# RadonShadow — Figure Framework Design & Text-Only Candidate Prompts

> **Paper:** Deterministic Quantum Shadow Tomography via MUB-DPRT Algebraic Duality
> **Visual Style:** Academic journal illustration, flat vector style, blue-primary color palette (`#1A5276` deep navy, `#2980B9` mid-blue, `#85C1E9` light blue, `#D4E6F1` accent), clean white background, thin black axis lines, IEEE double-column width (3.5 in × 2–3 in), sans-serif labels (Helvetica-style), no photorealistic elements.
> **Global Constraint:** All prompts produce **text-only schematic vector illustrations**. Mathematical notation rendered as clean LaTeX-style glyphs. No photorealistic 3D, no isometric, no skeuomorphic shadows.

---

## FIG1: Conceptual Overview — Breaking the Randomness Myth

**Panel Role:** First-glance reader hook. Instantly communicates that deterministic MUB measurement (=DPRT) is algebraically equivalent to random Pauli sampling but structurally superior.

**Layout:** Top–Bottom split with a bold central `=` connector.

**Top Half — "Random Sampling Paradigm (HKP)"**
- Left: A stylized **six-sided die** (flat icon, cube with dots) with a question mark on one face, symbolizing probabilistic choice.
- Center: A circular **MUB roulette wheel** divided into d+1 colored sectors (subtle blue gradient), with a pointer frozen mid-spin. Labels around the rim read `MUB₀`, `MUB₁`, …, `MUB_d` in small sans-serif.
- Right: A translucent **Pauli expectation bar chart** (3 columns: `⟨X⟩`, `⟨Y⟩`, `⟨Z⟩`) with wide error bars, annotated "`Var ∝ 1/M`".
- Caption ribbon across the top: *"Random Clifford Sampling: v ∝ 1/M (Shadow Bound)"*

**Middle — The Equality** (center, bold, spans full width)
- A large orange-red `≡` symbol inside a rounded rectangle, flanked by the labels *"Algebraic Duality"* left and *"T1 Theorem"* right.
- Thin arrows from the upper boxes converge into this `≡`, then descend.

**Bottom Half — "Deterministic DPRT Paradigm (RadonShadow)"**
- Left: A **discrete Radon grid** (4×4 lattice of dots with p+1 fixed lines overlaid, each line a different blue shade). The lines are perfectly regular — no randomness.
- Center: A **DPRT sinogram** (a small heatmap matrix: d rows, p+1 columns, blue-to-white color scale), annotated "d × (p+1) sinogram". The regularity of the pattern is visually striking compared to the chaotic upper panel.
- Right: A **variance comparison bar** — a single tall blue bar labeled "`Var ∝ 1/d²`" next to the faded Pauli bars from the upper panel, with a dashed arrow and annotation "**400–500× faster**".
- Caption ribbon across the bottom: *"Deterministic DPRT: Full p+1 directions, d²-fold variance reduction"*

**Design Notes:**
- The visual contrast between "messy top" (dice, roulette, scattered bars) and "ordered bottom" (grid, structured sinogram, clean bars) is the key rhetorical device.
- The orange `≡` symbol is the visual anchor — readers' eyes should land there.

### Prompt (English, Stable Diffusion / DALL·E 3)

```
Academic journal figure, flat vector illustration, IEEE double-column width 3.5 inches, clean white background, blue color scheme. Top half: a cartoon six-sided die next to a circular roulette wheel divided into d+1 colored sectors labeled MUB_0 through MUB_d, pointer mid-spin, with a Pauli expectation bar chart showing X, Y, Z columns with wide error bars. Bold caption reads "Random Clifford Sampling: Var ∝ 1/M". Center: a large orange-red triple-equals symbol inside a rounded box labeled "Algebraic Duality — T1 Theorem". Arrows converge from top panels into the equals sign. Bottom half: a 4x4 discrete grid of dots with p+1 neat diagonal lines overlaid in different blue shades (deterministic Radon lines), a small heatmap sinogram matrix labeled "d × (p+1)", and a side-by-side bar comparison showing a tall blue bar labeled "Var ∝ 1/d²" dwarfing the faded Pauli bars with annotation "400-500x faster". Bottom caption: "Deterministic DPRT: Full p+1 directions, d²-fold variance reduction". Sans-serif labels, no 3D, no shadows, clean academic vector style.
```

---

## FIG2: MUB-DPRT Algebraic Correspondence (T1 Theorem)

**Panel Role:** The mathematical heart of the paper. A three-column bridge diagram that proves MUB measurement statistics are algebraically identical to DPRT line integrals via a finite-field bijection.

**Layout:** Three vertical columns connected by labeled horizontal arrows.

**Column A — "MUB Measurement (Quantum)"**
- Top: Dirac ket notation box: `|ψ⟩ ∈ C^d` (state vector, light blue fill).
- Next: A copy-on-measure block: `{M_b^(j) = |u_b^(j)⟩⟨u_b^(j)|}_{j=0}^d` (MUB projectors for the j-th basis).
- Next: Probability formula block: `p_j(b) = ⟨ψ| M_b^(j) |ψ⟩ = |⟨u_b^(j) | ψ⟩|²` (Fourier sampling, d+1 bases × d outcomes).
- Bottom: Small grid icon showing d+1 rows of d probabilities → "`(d+1) × d` probability table."

**Column B — "Finite Field Bijection φ (p prime, d = p)"**
- Center: A **commutative diagram** (mini square):
  - Top-left: `F_p` (finite field)
  - Top-right: `{0,1,…,p−1} ⊕ {∞}` (extended index set, d+1 elements)
  - Bottom-left: `a ∈ F_p^*`  (nonzero elements = multiplicative group)
  - Bottom-right: `g^k mod p`  (primitive root powers, splitting into cosets)
- Diagonal arrows labeled: `φ: x ↦ log_g(x)` and `φ⁻¹`.
- Below the diagram: three small inline blocks:
  - `Weyl operators: X^a Z^b`
  - `Characters: χ_s(t) = ω_p^{st}`
  - `Trace map: Tr_F_p/F_1`
- Bottom annotation: *"Algebraic equivalence: MUB basis ↔ Radon projection direction"*

**Column C — "DPRT Sinogram (Classical)"**
- Top: A **discrete line integral formula**: `R_f(s, m) = Σ_{(i,j)∈L_{s,m}} f[i,j]` where `L_{s,m}` is the line with slope `s` and intercept `m`.
- Center: A rendered **sinogram heatmap** (d rows × (p+1) columns, blue-to-white), with a few example lines overlaid on a small 5×5 grid to the side showing how slope `s` and intercept `m` define a specific line.
- Bottom: Parameter counts: `p+1 slopes` (including ∞ for vertical), `p intercepts each`, total `(p+1)·p` measurements — exactly matching Column A's output.

**Arrows (horizontal, connecting columns):**
- Arrow A→B: labeled `"MUB projectors ↔ Weyl operators (via discrete Wigner)"`
- Arrow B→C: labeled `"Character sums ↔ Discrete line integrals (Fourier slice theorem over F_p)"`
- A broader bidirectional arrow spanning A↔C at the bottom: `"d+1 MUBs ≡ p+1 DPRT projections ≡ (d+1) × d probability table"`

### Prompt (English, Stable Diffusion / DALL·E 3)

```
Academic journal figure, flat vector illustration, IEEE double-column 3.5 inch width, clean white background, blue and navy color scheme, three-column layout. Left column titled "MUB Measurement (Quantum)": boxes showing Dirac ket notation |ψ⟩, MUB projector set {M_b^(j)}, probability formula p_j(b) = |⟨u_b^(j)|ψ⟩|², and a grid icon labeled "(d+1)×d probability table". Center column titled "Finite Field Bijection φ (F_p)": a commutative diagram square with nodes F_p, extended index set, F_p^*, and primitive root powers g^k mod p, connected by arrows labeled φ and φ⁻¹; below it show Weyl operators X^a Z^b, character sums χ_s(t), and trace map blocks. Right column titled "DPRT Sinogram (Discrete Radon)": formula R_f(s,m) = sum over line L_{s,m}, a sinogram heatmap matrix (d rows × p+1 columns, blue-to-white gradient), a small 5x5 grid with overlaid diagonal scan lines showing slope s and intercept m, and parameter count "(p+1) slopes × p intercepts". Horizontal arrows connecting columns: left-to-center labeled "MUB ↔ Weyl (discrete Wigner)", center-to-right labeled "Characters ↔ Line integrals (Fourier slice F_p)". Bottom bidirectional arrow spanning all three: "d+1 MUBs ≡ p+1 DPRT projections ≡ (d+1)×d probability table". All text in clean sans-serif, LaTeX-style math rendering, no 3D, no shadows, vector academic style.
```

---

## FIG3: Deterministic Shadow Protocol Pipeline

**Panel Role:** Algorithmic walkthrough. Shows the complete RadonShadow workflow from state preparation to observable estimation, with per-step complexity annotations. This is the "how it works" figure.

**Layout:** Linear pipeline — 5 connected stages arrayed horizontally (or 2-row if space demands: 3 top, 2 bottom). Each stage is a rounded rectangle with icon, label, and O(·) annotation.

**Stage 1 — "State Preparation"**
- Icon: Bloch sphere or a simple circle with `|ψ⟩` inside.
- Label: *"Prepare n copies of |ψ⟩"*
- Input arrow: "Input: ρ (unknown d-dim state)"
- Complexity badge: `O(1)` prep

**Stage 2 — "Deterministic Direction Selection"**
- Icon: Compass rose with exactly p+1 arrows radiating at fixed angles (not random — geometrically precise).
- Label: *"Select all p+1 MUB directions"*
- Sub-annotation: *"(HKP samples M directions; we fix M = p+1)"*
- Complexity badge: `O(p) = O(d)` selection

**Stage 3 — "Omnidirectional Projective Measurement"**
- Icon: A measurement apparatus (stylized detector) with p+1 light beams converging on a state.
- Label: *"Project onto |u_b^(j)⟩ for each direction j"*
- Sub-annotation: *"For each j: measure all d outcomes b ∈ {0,…,d−1}"*
- Complexity badge: `O(d²)` measurements

**Stage 4 — "DPRT Sinogram Construction"**
- Icon: A matrix being filled row by row (animated-still style — partially filled grid with data flowing in).
- Label: *"Build d × (p+1) sinogram S[j,b] = p_j(b)"*
- Sub-annotation: *"Each column = one projection direction; each row = one intercept"*
- Complexity badge: `O(d²)` assembly

**Stage 5 — "Fast Inversion & Estimation"**
- Icon: `S → inverse transform → {o_i}` — a sinogram with an arrow through a filter symbol (ramp filter icon) pointing to a list of observables.
- Label: *"DPRT inverse + median-of-means"*
- Sub-annotation: *"FBP: O(d² log d) or direct: O(d³)"*
- Output arrow: "Prediction: Tr(O_i ρ) ± ε"
- Complexity badge: `O(d² log d)` inversion

**Bottom summary bar:**
- A narrow strip across the full width: *"Total: O(d² log d) vs. HKP O(M d log M) with M ∝ log(d) — Deterministic amortization eliminates 1/√M sampling overhead."*

### Prompt (English, Stable Diffusion / DALL·E 3)

```
Academic journal figure, flat vector illustration, IEEE double-column 3.5 inch width, clean white background, blue-navy color scheme, horizontal pipeline of 5 connected stages. Stage 1 "State Preparation": Bloch sphere icon, label "Prepare n copies of |ψ⟩", O(1) badge. Stage 2 "Deterministic Direction Selection": compass rose icon with exactly p+1 fixed arrows at precise angles, label "Select all p+1 MUB directions", note "HKP samples M dirs; we fix M=p+1", O(d) badge. Stage 3 "Omnidirectional Projective Measurement": stylized detector icon with p+1 converging beams, label "Project onto |u_b^(j)⟩", note "all d outcomes per direction", O(d²) badge. Stage 4 "DPRT Sinogram Construction": partially filled matrix icon with data flowing in, label "Build d×(p+1) sinogram S[j,b]", note "columns=directions, rows=intercepts", O(d²) badge. Stage 5 "Fast Inversion & Estimation": sinogram icon with arrow through ramp filter symbol pointing to observable list, label "DPRT inverse + median-of-means", note "FBP: O(d² log d)", O(d² log d) badge. Arrows connect stages left-to-right. Bottom summary strip: "Total O(d² log d) vs HKP O(M d log M) — Deterministic amortization removes 1/√M overhead". All labels sans-serif, LaTeX math, no 3D, clean vector academic style.
```

---

## FIG4: Numerical Experiment Summary (4-Panel Grid)

**Panel Role:** Comprehensive experimental evidence in one figure. Four sub-panels in a 2×2 grid, each testing a different aspect of RadonShadow's performance.

**Layout:** 2×2 grid. Each sub-panel is labeled (a)–(d) at top-left. A summary statistics ribbon runs across the top.

**Top Summary Ribbon:**
- Three key numbers in large bold navy:
  - `median(λ) = 0.864`
  - `63.7% of dimensions < 0.90`
  - `400–500× speedup vs HKP`
- Small annotation: *"Across 168 prime dimensions d ∈ [3, 2039]"*

**Sub-panel (a) — "1-Local Observable Ratio λ vs. Dimension d"**
- Type: **Scatter plot** with log-scale x-axis.
- X-axis: "Dimension d (log scale)" from 3 to 2039.
- Y-axis: "1-Local Ratio λ = Var(HKP)/Var(RadonShadow)" from 0 to 1.2.
- Data: 168 blue dots (one per prime dimension).
- Reference line: black dashed horizontal at `λ = 1.0` labeled "parity".
- Reference line: orange dashed horizontal at `λ = 0.864` labeled "median".
- Trend curve: navy loess/LOWESS smooth fit through the scatter.
- Region shading: above `λ=1` is light blue ("RadonShadow Wins"), below `λ=1` is light gray ("Statistical Tie").
- Outlier annotation: point at d=2039 with a callout box: "d=2039: λ = 0.974, near-parity but still favorable".

**Sub-panel (b) — "Primitive Root Grouping Boxplot"**
- Type: **Grouped boxplot** (4 groups).
- X-axis: "Multiplicative Subgroup" with 4 tick labels: `H₁`, `H₂`, `H₄`, `H_full`.
- Y-axis: "1-Local Ratio λ".
- Each group: a box (IQR with median line) + whiskers + individual jittered dots behind.
- Color coding by subgroup size: darker blue = larger subgroup.
- Annotation: "Smallest subgroup H₁ → tightest variance" with a downward arrow from H₁ box.
- Small table below x-axis: `|H₁|=(p-1)/4, |H₂|=(p-1)/2, |H₄|=p-1` showing subgroup cardinalities.

**Sub-panel (c) — "Composite-Dimension Heatmap"**
- Type: **Heatmap** (square matrix). 
- X-axis: "Factor q₁" (composite factors), Y-axis: "Factor q₂".
- Cells: color intensity = `λ(d)` for `d = q₁·q₂` where q₁, q₂ are small primes/near-primes.
- Color scale: blue (λ ≪ 1, RadonShadow dominant) → white (λ ≈ 1, tie) → red (λ > 1, HKP wins).
- Annotation arrows pointing to blue regions: *"prime×prime: λ ≪ 1"*, to red regions: *"product of large co-prime factors"*.
- Diagonal: `q₁ = q₂` cells outlined in bold — annotated "pure prime powers: favorable".

**Sub-panel (d) — "Noise Robustness Curve"**
- Type: **Line plot** with error bands.
- X-axis: "Depolarizing Noise Strength η" from 0 to 1.
- Y-axis: "Relative Estimation Error ‖ô − o‖₂".
- Two curves:
  - Blue solid: "RadonShadow" with narrow shaded error band (±1 std).
  - Gray dashed: "HKP Shadow" with wider shaded error band (±1 std).
- Annotation: "RadonShadow maintains 2–3× lower error at all noise levels" with a brace spanning the gap between curves.
- Inset: small bar chart showing "Wall-clock time (ms)" comparison — RadonShadow bar much shorter than HKP bar.

### Prompt (English, Stable Diffusion / DALL·E 3)

```
Academic journal figure, flat vector illustration, IEEE double-column 3.5 inch width, clean white background, blue-navy color scheme, 2x2 grid layout with four sub-panels labeled (a)-(d). Top summary ribbon reads "median(λ)=0.864 | 63.7% dims < 0.90 | 400-500x speedup vs HKP" in bold navy.

Panel (a): Scatter plot, x-axis "Dimension d (log scale)" 3-2039, y-axis "1-Local Ratio λ" 0-1.2. 168 blue dots, gray dashed horizontal at λ=1.0 labeled "parity", orange dashed line at λ=0.864 "median", navy loess trend curve, light blue shading above λ=1 ("RadonShadow Wins"), light gray below. Outlier callout at d=2039: "λ=0.974".

Panel (b): Grouped boxplot, x-axis "Multiplicative Subgroup" with labels H₁, H₂, H₄, H_full, y-axis "λ". Boxes colored darker blue for smaller subgroups, individual jittered dots. Arrow: "Smallest H₁ → tightest variance". Table below: subgroup cardinalities.

Panel (c): Heatmap, x-axis "Factor q₁", y-axis "Factor q₂", color blue(λ≪1)→white(λ≈1)→red(λ>1). Arrows to blue regions: "prime×prime: λ≪1", to red regions: "large co-prime factors". Diagonal cells bold-outlined "pure prime powers: favorable".

Panel (d): Line plot, x-axis "Depolarizing Noise η" 0-1, y-axis "Relative Error ‖ô−o‖₂". Blue solid line "RadonShadow" with narrow error band, gray dashed "HKP Shadow" with wide error band. Brace annotation: "2-3x lower error at all η". Inset bar chart: wall-clock time comparison, RadonShadow bar far shorter.

All labels sans-serif, LaTeX math, no 3D, no shadows, clean vector academic style.
```

---

## FIG5: Prime Landscape — Dimension-Advantage Phase Diagram

**Panel Role:** The "hero shot" of the numerical results. Shows that across all 168 prime dimensions, RadonShadow dominates or ties HKP almost everywhere, with only sparse regions of HKP advantage.

**Layout:** Large central scatter plot with rich annotations.

**Main Plot:**
- X-axis: "Dimension d (= prime p)" from 0 to 2100, linear scale with tick marks at 500, 1000, 1500, 2000.
- Y-axis: "1-Local Ratio λ = Var(HKP)/Var(RadonShadow)" from 0.0 to 1.3.
- Data:
  - 168 **circle markers** (open circles, navy outline, light blue fill), one per prime dimension.
  - **Top-10 highlight**: 10 filled blue circles with larger radius at the dimensions with highest λ values (closest to HKP advantage). Each annotated with its `(d, λ)` value via a leader line.
- Reference lines:
  - `λ = 1.0`: thick black dashed line labeled "Break-even (λ = 1.0)".
  - `λ = 0.864`: thin orange dashed line labeled "Empirical median".
- **Trend fit**: A navy cubic-spline smooth curve through the scatter. Shows λ generally decreasing with d, asymptoting around 0.80–0.85.
- **Phase regions** (background shading):
  - Below `λ = 0.95`: **Light blue region** → labeled "DPRT Wins" with annotation: *"λ < 0.95: RadonShadow strictly dominates (≥5% advantage)"*
  - `λ ∈ [0.95, 1.05]`: **Light gray region** → labeled "Statistical Tie" with annotation: *"Indistinguishable given shot noise"*
  - Above `λ = 1.05`: **Very light pink region** → labeled "MUB Wins" (only ~3 primes fall here).
- **d=2039 outlier**: A red circle with a callout box: *"d = 2039: λ = 0.974. Largest prime tested. Near-parity — DPRT still slightly favorable."*
- **Small primes anomaly**: A dashed rectangle around d ∈ [3, 50] with annotation: *"Small-d regime: higher variance due to finite-size effects; DPRT advantage emerges for d > 50."*

**Side panel (right of main plot, narrow):**
- Small vertical bar chart: "Count of primes per phase region"
  - Blue bar: "DPRT Wins" — ~150 primes
  - Gray bar: "Statistical Tie" — ~15 primes
  - Pink bar: "MUB Wins" — ~3 primes

### Prompt (English, Stable Diffusion / DALL·E 3)

```
Academic journal figure, flat vector illustration, IEEE double-column 3.5 inch width, clean white background, blue-navy color scheme. Main scatter plot: x-axis "Dimension d (= prime p)" 0-2100 linear, y-axis "1-Local Ratio λ" 0.0-1.3. 168 open circle markers (navy outline, light blue fill), one per prime. 10 filled-blue larger circles highlighting Top-10 primes with leader-line callouts showing (d, λ) values. Thick black dashed line at λ=1.0 labeled "Break-even". Thin orange dashed line at λ=0.864 labeled "Empirical median". Navy cubic-spline trend curve showing λ decreasing asymptotically. Three phase regions as background shading: light blue below λ=0.95 ("DPRT Wins: λ<0.95, ≥5% advantage"), light gray λ∈[0.95,1.05] ("Statistical Tie"), very light pink above λ=1.05 ("MUB Wins"). Red circle at d=2039 with callout box: "d=2039: λ=0.974, largest prime, DPRT slightly favorable". Dashed rectangle around d∈[3,50] with annotation "Small-d: finite-size effects, DPRT advantage emerges d>50". Right side: narrow vertical bar chart "Count per Phase" with ~150 blue (DPRT Wins), ~15 gray (Tie), ~3 pink (MUB Wins). Sans-serif labels, LaTeX math, no 3D, clean vector academic style.
```

---

## FIG6: State Purity Effect — Variance Decomposition

**Panel Role:** Explains **why** deterministic measurement is especially powerful for pure states. Decomposes total variance into quantum (state) and classical (measurement) components, showing the former vanishes for pure states under deterministic protocols.

**Layout:** Side-by-side comparison with a connecting explanatory equation block.

**Left Panel — "Pure State |ψ⟩⟨ψ| (rank 1)"**
- Top icon: A clean Bloch sphere with a single surface point → "Tr(ρ²) = 1".
- **Variance decomposition bar (horizontal stacked bar, full width):**
  - Segment 1 (dark blue, ~10% of bar): labeled "Classical (measurement shot noise)" with value `σ²_shot`.
  - Segment 2 (light blue, ~90% of bar): labeled "Quantum (state projection)" with value `σ²_proj ≈ 0*`.
  - Asterisk footnote: "*Vanishes under deterministic full-basis projection"
- Below the bar: a formula block:
  ```
  Var(ô) = Var_shot(ô) + Var_proj(ô)  
         = (1/n)·Var_classical + 0
         → Scales as 1/n, not 1/(nM)
  ```
- Small annotation: *"For pure states, deterministic full MUB measurement eliminates projection variance entirely — the dominant noise source in random sampling."*

**Right Panel — "Mixed State ρ (rank > 1)"**
- Top icon: A fuzzy Bloch sphere (interior point, not surface) → "Tr(ρ²) < 1".
- **Variance decomposition bar (horizontal stacked bar, full width):**
  - Segment 1 (dark blue, ~40%): labeled "Classical (measurement shot noise)".
  - Segment 2 (medium blue, ~60%): labeled "Quantum (state projection, non-zero)" with value `σ²_proj > 0`.
- Below the bar: a formula block:
  ```
  Var(ô) = Var_shot(ô) + Var_proj(ô)  
         = (1/n)·Var_classical + (1/n)·f(ρ)
         → DPRT still better, but gap narrows with mixedness
  ```
- Small annotation: *"For mixed states, projection variance does not vanish; DPRT advantage persists but partially offset by state mixedness."*

**Center Connector (between the two panels):**
- A vertical strip containing:
  - A **purity axis** (arrow from top to bottom): labeled "State Purity Tr(ρ²)" with "1.0 (Pure)" at top and "1/d (Maximally Mixed)" at bottom.
  - A **ratio curve**: A simple line plot showing `Var(RadonShadow)/Var(HKP)` vs. purity — starting very low (<< 1) at purity=1, curving upward and approaching ~1 at purity→1/d.
  - Annotation on the curve: *"DPRT advantage ∝ purity. Pure states → maximal gain."*

### Prompt (English, Stable Diffusion / DALL·E 3)

```
Academic journal figure, flat vector illustration, IEEE double-column 3.5 inch width, clean white background, blue-navy color scheme, left-right split layout with center connector. Left panel titled "Pure State |ψ⟩⟨ψ| (rank 1)": Bloch sphere icon with surface point, horizontal stacked variance bar: dark blue "Classical σ²_shot ~10%", light blue "Quantum σ²_proj ≈ 0*" ~90%, asterisk note "vanishes under deterministic full-basis projection". Formula block: Var(ô) = Var_shot + Var_proj = (1/n)·Var_classical + 0 → scales as 1/n. Annotation: "Deterministic full MUB eliminates projection variance — the dominant noise source in random sampling." Right panel titled "Mixed State ρ (rank>1)": fuzzy Bloch sphere interior point, horizontal stacked variance bar: dark blue "Classical ~40%", medium blue "Quantum σ²_proj>0 ~60%". Formula block: Var(ô) = (1/n)·Var_classical + (1/n)·f(ρ) → DPRT still better but gap narrows. Annotation: "Projection variance nonzero; DPRT advantage partially offset by mixedness." Center vertical strip: purity axis arrow top-to-bottom "Tr(ρ²): 1.0 (Pure) → 1/d (Maximally Mixed)", with line plot showing Var(RadonShadow)/Var(HKP) ratio vs purity, starting <<1 at purity=1, curving to ~1 at purity→1/d. Annotation on curve: "DPRT advantage ∝ purity. Pure states → maximal gain." Sans-serif labels, LaTeX math, no 3D, clean vector academic style.
```

---

## Style Reference Card (Apply to All Figures)

```
GLOBAL STYLE DIRECTIVES:
- Color palette:
  Primary navy:    #1A5276 (deep blue, axes/labels/main elements)
  Mid blue:        #2980B9 (data markers, boxes, emphasis)
  Light blue:      #85C1E9 (fills, backgrounds, secondary elements)
  Accent blue:     #D4E6F1 (lightest fills, phase regions)
  Orange accent:   #E67E22 (callouts, equality symbol, highlights — use sparingly)
  Red accent:      #C0392B (outlier markers only)
  Gray:            #7F8C8D (reference lines, secondary text)
  Light gray:      #ECF0F1 (phase regions, neutral fills)
  Light pink:      #FADBD8 (MUB Wins phase region only)
- Typography: Helvetica or Arial equivalent, 8pt labels, 10pt titles, 7pt annotations
- Axes: Thin black (0.5pt), tick marks inward, no top/right spines where applicable
- Data markers: 1.5pt stroke, 3pt radius for circles
- Line weights: 1pt for curves, 0.75pt for reference lines, 2pt for emphasis lines
- All mathematical notation rendered as LaTeX (italic for variables, upright for operators)
- No gradients (use flat fills); no drop shadows; no 3D extrusion
- Clean, minimal, suitable for IEEE Transactions on Information Theory
```

---

## Figure Dependency Graph

```
FIG1 (Concept) ──────────────────────┐
                                     ├── FIG3 (Pipeline) ── FIG4 (Experiments)
FIG2 (Theorem) ──────────────────────┘        │
                                              ├── FIG5 (Prime Landscape)
                                              └── FIG6 (Purity Effect)
```

- **FIG1 + FIG2** are the foundational pair: intuition + mathematics. Place early (Sections 1–2).
- **FIG3** explains the protocol. Place in Section 3 (Methods).
- **FIG4** is the comprehensive results summary. Place in Section 4 (Results).
- **FIG5** is the detailed prime-dimension analysis. Place in Section 4.2 or Appendix.
- **FIG6** is the mechanistic explanation. Place in Section 5 (Discussion).

---

*Generated for RadonShadow v1.0 — Figure Design Framework*
*Date: 2026-07-15*
