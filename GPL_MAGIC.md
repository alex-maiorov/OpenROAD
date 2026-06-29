# GPL Magic Numbers

This document catalogs every hardcoded (non-configurable) constant, threshold,
and heuristic formula in `src/gpl/` that has a **real impact** on placement
quality, convergence behaviour, or runtime.  Configurable `PlaceOptions` 
defaults are noted only when they interact with hardcoded logic; purely 
Tcl-configurable knobs (e.g. `density`, `overflow`) are omitted.

---

## 1. NesterovPlace — Core Loop & Coefficients

### 1.1 Wire-length coefficient (`updateWireLengthCoef`)
**File:** `src/gpl/src/nesterovPlace.cpp:1340–1354`
**What it does:** Maps the normalised overflow (0–∞) onto a wire-length vs.
density trade-off coefficient.
**Hardcoded formula:**
| Condition       | `wireLengthCoef` |
|-----------------|-------------------|
| `overflow > 1`  | `0.1`             |
| `overflow < 0.1`| `10.0`            |
| otherwise       | `1.0 / 10^((overflow - 0.1) * 20/9 - 1)` |
**Impact:** Directly determines how strongly the placer pulls instances toward 
wire-length-minimising positions vs. spreading them out.  The numbers `0.1`, 
`10.0`, `20`, `9`, and the exponent base `10.0` are all hardcoded and 
cannot be tuned from Tcl.
**Final multiplier:** `wireLengthCoefX_ *= baseWireLengthCoef_` (and 
analogous for Y).

### 1.2 Nesterov acceleration coefficient
**File:** `src/gpl/src/nesterovPlace.cpp:1184`
```
curA = (1.0 + sqrt(4.0 * prevA * prevA + 1.0)) * 0.5
```
**Impact:** This is the textbook Nesterov momentum update (Algorithm 1 in 
ePlace-MS).  The constants `1.0`, `4.0`, `0.5` are the canonical formula.

### 1.3 Recursion limits for zero-gradient
**File:** `src/gpl/src/nesterovBase.h:821–822` (`NesterovPlaceVars`)
```
static constexpr int maxRecursionWlCoef   = 10;
static constexpr int maxRecursionInitSLPCoef = 10;
```
**Impact:** When wire-length gradient sum is exactly zero (e.g., all pins 
overlap), the code halves `wireLengthCoefX_`, `wireLengthCoefY_`, and 
`baseWireLengthCoef_` and retries up to 10 times.  When the 
initial step-length is NaN/Inf, `initialPrevCoordiUpdateCoef` is multiplied 
by 10 and retried.  Both are hard limits preventing infinite recursion.

### 1.4 Backtracking
**File:** `src/gpl/src/nesterovBase.h:811` (`NesterovPlaceVars`)
```
static constexpr int maxBackTrack = 10;          // nesterovPlaceVars
```
**File:** `src/gpl/src/nesterovPlace.cpp:826–880`  (`doBackTracking`)
**Impact:** During each Nesterov iteration the placer tries up to 10 backtrack 
steps.  On each step it recomputes the gradient and step-length.  Backtracking
continues while the new step-length shrinks to ≤ 95 % of the current
(and stays ≥ 0.01); it stops once the step-length fails to shrink farther,
bottoms out below 0.01, or becomes NaN/Inf.
**Hardcoded threshold** inside `nesterovUpdateStepLength` 
(`nesterovBase.cpp:3784–3801`):
- `newStepLength > stepLength_ * 0.95` → accept (exit backtrack)
- `newStepLength < 0.01` → floor at `0.01`
- `nan` / `inf` → divergent

### 1.5 Minimum preconditioner
**File:** `src/gpl/src/nesterovBase.h:815`
```
static constexpr float minPreconditioner = 1.0;
```
**Impact:** Applied in **two** gradient-update functions to ensure the summed 
preconditioner never drops below 1, preventing gradient blow-up:

- `updateGradients` (bulk, `nesterovBase.cpp:3346–3349`)
- `updateSingleGradient` (per-cell backtracking, `nesterovBase.cpp:3538–3539`)

Both locations use the same logic:
```cpp
sumPrecondi.x = std::max(sumPrecondi.x, NesterovPlaceVars::minPreconditioner);
```

### 1.6 Initial previous-coordinate update coefficient
**File:** `src/gpl/src/nesterovBase.h:816`
```
float initialPrevCoordiUpdateCoef = 100;
```
**Impact:** Sets the distance between `prevSLPCoordi` and `curSLPCoordi` 
during `init()` (`updateInitialPrevSLPCoordi`, `nesterovBase.cpp:3545–3570`).
Multiplied by 10 on retry when the initial step-length is NaN/Inf 
(`nesterovPlace.cpp:252`).  This value is hardcoded in `NesterovPlaceVars` 
and is **not** configurable from Tcl (`PlaceOptions` has no corresponding 
field, and the `NesterovPlaceVars` constructor does not read it from 
options).

### 1.7 Divergence snapshot overflow threshold
**File:** `src/gpl/src/nesterovPlace.cpp:1377`
```cpp
if (hpwl < min_hpwl_ && average_overflow_unscaled_ <= 0.25)
```
**Impact:** A divergence-revert snapshot is only saved when overflow is 
≤ 0.25.  This protects against saving a snapshot during early high-overflow 
iterations.

### 1.8 GIF/Image debug hardcodes
**File:** `src/gpl/src/nesterovPlace.cpp` (multiple locations)
| Constant | Value | Where |
|----------|-------|-------|
| GIF width (px) | `500` | lines 338, 612, 734 |
| dbu_per_pixel divisor | `1000.0` | lines 341, 615, 738 |
| GIF frame delay | `20` | lines 342, 616, 739 |
| GIF save stride | `iter % 10 == 0` | line 336 |

---

## 2. NesterovBase — Density & Bin Grid

### 2.1 Density size floor / sqrt(2) rule
**File:** `src/gpl/src/nesterovBase.cpp:46`
```cpp
#define REPLACE_SQRT2 1.414213562373095048801L
```
**Impact:** In `updateDensitySize` (lines 2547–2576): if a GCell's physical 
dimension is smaller than `sqrt(2) * binSize`, its *density* size is bumped 
up to `sqrt(2) * binSize` and its `densityScale` is reduced proportionally.
This prevents cells from being invisible to the density grid.

### 2.2 Filler-cell statistics window
**File:** `src/gpl/src/nesterovBase.cpp:2156–2157`
```cpp
int minIdx = dxStor.size() * 0.05;
int maxIdx = dxStor.size() * 0.95;
```
**Impact:** Filler cells (virtual cells that fill white-space for density 
modelling) are sized to the average of the 5 %–95 % percentile range of 
real instance dimensions.  This discards outliers.

### 2.3 Filler-count limit
**File:** `src/gpl/src/nesterovBase.cpp:2230–2268`
```cpp
const double limit_filler_ratio = 10;
const double max_edge_fillers = 1024;
```
**Impact:** If the number of fillers exceeds 10× the number of real instances, 
the filler cell size is scaled up.  After scaling, filler dimensions are
capped at `max(core_size / 1024, original_filler_size)`:
if the original filler was smaller than `core/1024`, it may at most grow to
`core/1024`; if already larger, it may not grow at all.  This keeps at
least ~1024 fillers across the core edge and prevents fillers that are
already large from becoming even larger.

### 2.4 Random seed for initial perturbation
**File:** `src/gpl/src/nesterovBase.cpp:2040`
```cpp
srand(42);
```
**Impact:** All region NesterovBases use a fixed initial random seed.  Instance 
positions get a random offset `rand() % (2 * dbu_per_micron) - dbu_per_micron` 
(lines 2054–2055).  Because `rand()` is used (not `std::mt19937` like for 
fillers), the perturbation is deterministic but limited by `RAND_MAX=32767`.

### 2.5 Macro Gaussian smoothing
**File:** `src/gpl/src/nesterovBase.cpp:5042–5045`
```cpp
const biNormalParameters i = {meanX, meanY,
                               meanX / 6,   // sigmaX
                               meanY / 6,   // sigmaY
                               ...};
```
**Impact:** The overlap area between a macro and a bin is smoothed with a 
bivariate normal CDF whose σ = mean/6.  The divisor `6` controls how 
"smeared" the macro's blockage appears.

**File:** `src/gpl/src/nesterovBase.cpp:5061–5063`  (same function)
```cpp
if (scaled >= original) {
    return std::min<float>(scaled, original * 1.10);
}
```
**Impact:** Within `getOverlapArea`, the per-instance bivariate-normal smoothed 
overlap is capped at `1.10 ×` the geometric overlap.  The smoothed value is 
never allowed to exceed 110 % of the geometric rectangle; if it falls below 
the geometric overlap, the geometric value is returned instead.  This is a 
separate, per-macro cap from the bin-level aggregation cap in § 2.6.

### 2.6 Non-place-area cap multiplier
**File:** `src/gpl/src/nesterovBase.cpp:928`
```cpp
const int64_t cap = unionArea[i] * bin.getTargetDensity() * 1.10f;
```
**Impact:** When multiple fixed macros overlap the same bin, their combined 
blockage is capped at `1.10 × targetDensity ×` the geometric union area to 
prevent over-counting.  (The `targetDensity` factor — typically 0.6–0.9 — 
means the cap is tighter than `1.10 × unionArea` alone would imply.)

### 2.7 Overflow denominator guard
**File:** `src/gpl/src/nesterovBase.cpp:3656–3660`
```cpp
const float overflowDenominator =
    std::max(static_cast<float>(getNesterovInstsArea()),
             fractionOfMaxIters * pb_->nonPlaceInstsArea() * 0.05f);
```
**Impact:** In macro-dominated designs the denominator for overflow can be 
tiny, preventing convergence.  The denominator is floored at 
`(iter/maxIter) * nonPlaceInstsArea * 0.05` to create a gradual ramp.

### 2.8 Phi coefficient adjustment
**File:** `src/gpl/src/nesterovBase.cpp:3845–3862`  (`nesterovAdjustPhi`)
```cpp
if (!nbVars_.isMaxPhiCoefChanged && sum_overflow_unscaled_ < 0.35f) {
    nbVars_.isMaxPhiCoefChanged = true;
    maxPhiCoef *= 0.99;
}
if (maxPhiCoef <= 1.0f) {
    maxPhiCoef = 1.01f;
}
```
**Impact:** The **first time** overflow drops below `0.35`, `maxPhiCoef` is 
shrunk by factor `0.99` to improve convergence in large designs.  The 
`isMaxPhiCoefChanged` guard ensures this is a one-time adjustment — 
subsequent iterations will not shrink it further.  `maxPhiCoef` is never 
allowed below `1.01` (density penalty can only grow, never shrink).

### 2.9 minSumOverflow tracking
**File:** `src/gpl/src/nesterovBase.cpp:3596–3597`
```cpp
minSumOverflow_ = 1e30;
hpwlWithMinSumOverflow_ = 1e30;
```
**File:** `src/gpl/src/nesterovBase.cpp:3767`
```cpp
if (iter > 50 && minSumOverflow_ > sum_overflow_unscaled_) { ... }
```
**Impact:** The "minimum observed overflow" tracking only starts after 
iteration `50` to avoid recording the early high-overflow phase.

### 2.10 Wire-length preconditioner
**File:** `src/gpl/src/nesterovBase.cpp:1599–1603`
```cpp
NesterovBaseCommon::getWireLengthPreconditioner(const GCell*) const
{
  return FloatPoint(1.0f, 1.0f);
}
```
**Impact:** The wire-length preconditioner is hardcoded to `(1, 1)`.  This 
means the WL gradient is never scaled down based on cell size — only the 
density, timing, and routability preconditioners contribute to the 
denominator.

### 2.11 Weighted-average wire-length model floor
**File:** `src/gpl/src/nesterovBase.h:803`
```cpp
static constexpr float minWireLengthForceBar = -300;
```
**Impact:** In `updateWireLengthForceWA`, pin contributions where 
`exp(x/gamma)` is below `exp(-300)` are silently skipped.  This prevents 
underflow and speeds up the WA calculation for distant pins.

---

## 3. Divergence Detection

### 3.1 Primary divergence check
**File:** `src/gpl/src/nesterovBase.cpp:3995–4038`  (`checkDivergence`)
| Condition | Hardcoded Value |
|-----------|-----------------|
| Overflow must be below | `0.2` |
| Overflow absolute increase | `≥ 0.02` |
| HPWL increase factor | `≥ 1.2×` |
**Impact:** Divergence is flagged only when overflow is already low 
(`< 0.2`), has increased by at least 0.02, and HPWL has grown by 20 % 
relative to the min-observed point.

### 3.2 Secondary divergence check (reported values)
| Condition | Hardcoded Value |
|-----------|-----------------|
| Minimum overflow precondition | `minSumOverflow_ < 0.2` |
| Overflow change tolerance | `0.05` |
| HPWL increase tolerance | `0.25` (25 %) |
**Impact:** When `minSumOverflow_` is already below `0.2`, if overflow and
HPWL both increase by these thresholds compared to the 10-iteration rolling
reference, divergence is flagged.  This provides a second, coarser detection
path that does not rely on the per-iteration min-overflow snapshot point.

### 3.3 Timing-aware divergence suppression
**File:** `src/gpl/src/nesterovBase.cpp:4001–4002, 4007, 4026`  (and `:5815–5835`)
**Impact:** In timing-driven mode, both the primary (3.1) and secondary (3.2)
divergence checks are **additionally gated**: divergence is only flagged when
**timing has also degraded**, not solely on HPWL/overflow change.  This
prevents false-positive divergence exits when the solver legitimately trades
wirelength for timing improvements.

The gate is activated when `timingDrivenMode` is true and STA has been
queried at least once (`sta_update_count_ > 0`).  In that case, the
`timing_degraded_` flag must be true for divergence to be raised (lines
4007, 4026).  `timing_degraded_` is set to true in `updateSTA()` (line
5831) when **both** WNS and TNS degrade (become more negative) between
consecutive STA queries:
```cpp
timing_degraded_ = (current_wns < prev_wns_) && (current_tns < prev_tns_);
```
A single metric improving is enough to treat the iteration as a legitimate
timing-for-WL trade-off and suppress the divergence signal.

The state variables (`prev_wns_`, `prev_tns_`, `timing_degraded_`,
`sta_update_count_`) are tracked per-region in `NesterovBase`
(`nesterovBase.h:1421–1424`).

### 3.4 Consecutive divergence gating
**File:** `src/gpl/src/nesterovPlace.cpp:1232–1279`
When `isDiverged` returns true but `consecutiveDivergeCount_` hasn't reached
`divergeConsecutiveThreshold` (configurable, default `1`), the divergence 
state is reset and the solver continues.  This single-iteration tolerance 
prevents spurious exits.

---

## 4. RouteBase — Routability

### 4.1 Routability iteration limits
**File:** `src/gpl/src/routeBase.h:206–207`
```cpp
int max_routability_no_improvement_ = 3;
int max_routability_revert_ = 50;
```
**Impact:** Routability reverts stop after 3 consecutive non-improving 
iterations or 50 total reverts.  These are hardcoded and not configurable 
from Tcl.

### 4.2 Minimum congestion improvement threshold
**File:** `src/gpl/src/routeBase.cpp:578`
```cpp
if ((minRc_ - curRc) > 0.001)
```
**Impact:** An improvement in routing congestion must exceed `0.001` to be 
recognised as "better than minimum".

### 4.3 Congestion averaging percentiles
**File:** `src/gpl/src/routeBase.cpp:849–861`
```cpp
0.005, 0.01, 0.02, 0.05  // top 0.5%, 1%, 2%, 5%
```
**Impact:** The weighted routing congestion metric is computed from the 
average of the top 0.5 %, 1 %, 2 %, and 5 % most-congested tiles.  
The percentages are hardcoded; only the weighting coefficients (rcK1–rcK4) 
are configurable.

### 4.4 Tile blockage ignore ratio
**File:** `src/gpl/src/routeBase.cpp:225`
```cpp
float ignoreEdgeRatio = 0.8;
```
**Impact:** Tiles whose blockage exceeds 80 % of capacity are excluded from 
congestion calculations.

### 4.5 Minimum inflation ratio trigger
**File:** `src/gpl/src/routeBase.cpp:226`
```cpp
float minInflationRatio = 1.01;
```
**Impact:** Inflation is only applied to tiles whose congestion ratio ≥ 1.01.

### 4.6 RUDY congestion clamp
**File:** `src/gpl/src/nesterovBase.cpp:2794`
```cpp
std::isfinite(ratio) ? std::clamp(ratio, 0.0f, 10.0f) : 0.0f
```
**Impact:** Finite RUDY congestion ratios are clamped to `[0, 10]` before being 
stored for the routability gradient pass.  Non-finite (NaN/Inf) ratios are
replaced with `0.0f`.

---

## 5. Timing Gradient Pass

### 5.1 Slack query lower bound
**File:** `src/gpl/src/nesterovBase.cpp:5164`
```cpp
float slack_min = -1e30f;  // essentially -inf
```
**Impact:** When querying STA for violating paths, no lower bound is 
applied — all paths with slack ≤ `slack_upper` are captured.

### 5.2 Minimum slack threshold
**File:** `src/gpl/src/nesterovBase.h:1479`
```cpp
static constexpr float kMinSlackThreshold_ = 1e-3f;
```
**Impact:** Paths whose slack magnitude exceeds this threshold are skipped 
for per-cell timing gradient calculation (`std::abs(path.slack) > 1e-3` → 
`continue`).  The check is active at `nesterovBase.cpp:5596` and `:5718`;
commented out at `:5429` for the all-to-all path.

**Practical note:**  STA slack is in seconds, so the threshold value 
`1e-3` = 1 ms.  Real timing slacks are typically in the picosecond range 
(1e-10 to 1e-11 s), meaning `abs(slack) > 1e-3` is **never true** in 
practice — this filter is effectively a no-op for any real design.

### 5.3 λ₅ (L₅) timing-to-wirelength warning
**File:** `src/gpl/src/nesterovBase.cpp:3297`
```cpp
const float timing_to_wirelength_warning_thresh = 100.0;
```
**Impact:** If the mean timing gradient magnitude is ≥ 100× the mean 
wire-length gradient magnitude, a warning is emitted (code uses `>=` at 
line 3386).

---

## 6. InitialPlace — Conjugate Gradient

### 6.1 Solver convergence
**File:** `src/gpl/src/initialPlace.cpp:114`
```cpp
if (error_max <= 1e-5 && iter >= 5)
```
**Impact:** The CG solver stops when residual ≤ `1e-5` **and** at least 5 
iterations have completed.  Both thresholds are hardcoded.

### 6.2 Trigram matrix reserve
**File:** `src/gpl/src/initialPlace.cpp:275–276`
```cpp
listX.reserve(1'000'000);
listY.reserve(1'000'000);
```
**Impact:** Pre-allocation for the sparse matrix construction — may cause 
re-allocations on very large designs.

---

## 7. Miscellaneous

### 7.1 Macro detection row threshold
**File:** `src/gpl/src/placerBase.cpp:68`
```cpp
constexpr int row_limit = 6;
```
**Impact:** Instances taller than 6 standard-cell rows are classified as 
macros (if not already marked as BLOCK in LEF).

### 7.2 Bin count search
**File:** `src/gpl/src/nesterovBase.cpp:772–790`
- Starting point: `foundBinCnt = 2`
- Upper bound: `foundBinCnt <= 1024`, growing by factor `2`
- Minimum ideal bin count: `idealBinCnt = std::max(idealBinCnt, 4)`
**Impact:** The bin grid dimensions are powers of 2 between 2 and 1024, 
adjusted for aspect ratio.  The `1024` cap and `4` minimum constrain bin 
counts.

### 7.3 fastExp approximation
**File:** `src/gpl/src/nesterovBase.cpp:5097–5111`
```cpp
exp = 1.0f + exp / 1024.0f;
exp *= exp;  // 10 times
```
**Impact:** A 10th-order Taylor-inspired approximation of `exp(x)`.  The 
`1024` divisor and `10` iterations trade accuracy for speed in the WA 
wire-length model.

### 7.4 Timing gradient blend default
**File:** `src/gpl/include/gpl/Replace.h:113`
```cpp
float timingGradPassBlend = 0.3F;
```
**Impact:** `30 %` new gradient + `70 %` previous gradient blending 
prevents discontinuities when STA re-query changes the violating-path set.

### 7.5 Routability gradient default range cross-over
**File:** `src/gpl/include/gpl/Replace.h:129–130`
```cpp
int routabilityGradPassFirstIter = 800;
int routabilityGradPassRunInterval = 100;
```
**Impact:** The routability gradient pass starts at iteration 800 and 
refreshes congestion every 100 iterations by default.  These are 
configurable but noted because changing them dramatically alters the 
routability gradient's influence.

### 7.6 Path saturation parameters
**File:** `src/gpl/src/nesterovBase.cpp:5517–5521`
```cpp
const float L = std::max(path_span / nbVars_.timing_pass_saturation_kL,
                          nbVars_.timing_pass_saturation_minL);
```
**Impact:** Pseudo-Huber saturation characteristic length `L` uses 
`path_span / kL` (default kL=3.0) floored at `minL` (default 1000.0 DBU).
Short paths see nearly-linear forces; long paths saturate.

### 7.7 Instance overlap in checkConsistency
**File:** `src/gpl/src/nesterovBase.cpp:2388`
```cpp
const int64_t tolerance = 10000;
```
**Impact:** Consistency checks use a hardcoded area tolerance of 10,000 DBU².

### 7.8 Uniform density rounding
**File:** `src/gpl/src/nesterovBase.cpp:2215`
```cpp
uniformTargetDensity_ = ceilf(uniformTargetDensity_ * 100) / 100;
```
**Impact:** Uniform density is rounded up to two decimal places.

### 7.9 Incremental placement defaults
**File:** `src/gpl/src/replace.cpp:177–197`
```cpp
locked_options.overflow = std::max(options.overflow, 0.2f);
locked_options.nesterovPlaceMaxIter = 300;
final_options.initDensityPenaltyFactor = 1;
```
**Impact:** In incremental mode, overflow is floored at `0.2`, max iterations 
capped at `300`, and the final tightening pass uses `initDensityPenalty = 1`.

### 7.10 Site size defaults (PlacerBase)
**File:** `src/gpl/src/placerBase.h:416–417`
```cpp
int siteSizeX_ = 0;
int siteSizeY_ = 9;
```
**Impact:** Default site height of 9 DBU (overridden during init from DB).

### 7.11 Uniform target density epsilon
**File:** `src/gpl/src/nesterovBase.cpp:2196–2199`
```cpp
float tmp_targetDensity
    = static_cast<float>(stdInstsArea_)
          / static_cast<float>(whiteSpaceArea_ - macroInstsArea_)
      + 0.01;
```
**Impact:** When `useUniformTargetDensity` is enabled (which is the default for
incremental placement), the computed uniform density is bumped by `0.01`
(1 percentage point) to ensure the placer always has a small amount of extra
white space to work with, preventing convergence stalls at exactly 100 %
utilisation.

### 7.12 Pin-density instance scaling caps
**File:** `src/gpl/src/placerBase.cpp:897–900`
```cpp
if (scale > 1.2) {
    scale = 1.2;
} else if (scale < 0.95) {
    scale = 0.95;
}
```
**Impact:** When instances are resized to match average pin density (disabled
via `disablePinDensityAdjust`), the scale factor is clamped to `[0.95, 1.20]`.
This prevents any single instance from shrinking below 95 % or growing beyond
120 % of its original area, limiting downstream routability inflation effects.

