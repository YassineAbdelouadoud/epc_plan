# Linearized additive calibration model — assessment

*Context: exploration of replacing the multiplicative main-effects model of the EPC calibration with (a) a linearized model with additive factors, and (b) a two-component formulation in which the non-EPC energy (base load) is fitted alongside as an additive model with its own predictors.*

Both halves of this idea are good, and they're good in different ways: the linearization is a major tractability and *rigor* win at (almost) no interpretive cost, while the two-component idea — EPC-corrected consumption plus a separately-modeled non-EPC base load — is arguably a fix to a real structural misspecification in the current channels, not just a convenience.

## 1. The linearized model

Replace the product with a sum of deviations:

$$
\hat{c}_d \;=\; \Big(g_0 + \sum_v \beta_v[\ell_v(d)]\Big)\, E_d,
\qquad
\text{mass-weighted } \sum_k \beta_v[k] = 0 \ \text{ per variable } v.
$$

The predicted district total becomes **linear in the parameters**:

$$
\hat y_i \;=\; g_0\, E_i \;+\; \sum_v \sum_k \beta_v[k]\, M_{i,v,k},
$$

where

$$
M_{i,v,k} \;=\; (A S_v)_{i,k}
$$

is the conventional-consumption mass of level $k$ of variable $v$ in district $i$ — quantities already computed in the pipeline.

The parameter count is unchanged (one coefficient per level), and each coefficient still reads as a single statement, just in percentage-point rather than percent terms: "owner-occupied adds $+0.10$ to the correction factor" instead of "$\times\,1.19$". Since the fitted multipliers mostly live in $[0.8,\ 1.2]$, the two families differ only at second order,

$$
\prod_v (1+\delta_v) \;\approx\; 1 + \sum_v \delta_v \;+\; \underbrace{\sum_{v<w}\delta_v \delta_w}_{\text{a few percent}},
$$

so near-identical district-level fit is expected — cheap to verify with the existing CV harness.

### What tractability actually buys — more than speed

- **The design collapses.** The ~19,000-cell matrix is never needed for fitting, only the ~50–60 district × level mass columns. The fit is a weighted least squares on a 40,000 × ~60 matrix — milliseconds, with a closed form under ridge.
- **The BCD caveat paragraph disappears.** The methods section currently has to admit the block-coordinate elastic net "minimizes no single joint objective" and claims no convergence guarantee, with grouping properties holding only under an unmet orthonormality assumption. A linear model makes the group elastic net *exact*: one jointly convex problem, global optimum, and Zou & Hastie's grouping result applies literally rather than per-block-approximately. The renormalize-after-each-block heuristic is replaced by an exact linear constraint (mass-weighted zero mean per variable), which also resolves the dummy trap since $\sum_k M_{i,v,k} = E_i$ for every $v$.
- **Everything downstream gets cheap and standard.** The $\lambda_1$ path, leave-departments-out CV, the 200-replicate department bootstrap, and the per-candidate-split refits of the MOB extension all become trivial refits. MOB itself becomes the *textbook* case — linear regression in the leaves, exactly Zeileis et al.'s primary illustration, with standard scores instead of hand-derived nonlinear ones.
- **Fit and identification bounds live in one world.** The identification-bounds SOCP machinery is unchanged (it is model-agnostic over the free-cell class), but the fitted model is now itself a linear point in the same geometry, and cluster-robust standard errors come in closed form.

### The one thing the multiplicative form does that additive cannot: the occupancy gate

"Vacant $\Rightarrow$ consumption $\approx 0$ *whatever else the dwelling is*" is a multiply-by-zero statement; no additive $\beta_{\text{vacant}}$ can cancel $g_0$ plus deviations that vary cell by cell. Related: non-negativity of $\hat c_d$ is no longer automatic — a stack of negative deviations can push a cell's factor below zero. Both have clean fixes, best combined:

1. **Fixed multiplicative gates for vacant/occasional**, at their prior values — honest, since the identification analysis shows these levels are prior-anchored, not data-estimated. The gate just pre-multiplies $E_d$ before the linear fit.
2. **Cell-level non-negativity** $N\theta \ge 0$ **by delayed constraint generation** — solve unconstrained, check the ~19k observed cells, add constraints only if violated (likely rarely, once the gates handle occupancy). The typical solve stays closed-form.

## 2. The two-component model

This is the stronger half of the idea. Currently the non-EPC energy (cooking-gas floor; the non-heating electricity floor) enters as a **national-average constant folded into $E_d$**, which has two consequences: its spatial variation is ignored, and — worse — the fitted multipliers *scale it*. The tenure multiplier of 1.22 on the electric channel scales appliances and lighting exactly as much as heating, which is physically wrong: prebound is a heating phenomenon. On the electric channel this is not a rounding issue — the non-heating floor is roughly half the mass.

The proposal fixes this structurally, and it stays one linear model:

$$
\hat y_i \;=\;
\underbrace{\;g_0\, E_i + \sum_v \sum_k \beta_v[k]\, M_{i,v,k}\;}_{\substack{\text{EPC heating correction} \\ \text{(MWh-mass features)}}}
\;+\;
\underbrace{\;\alpha_0\, n_i + \sum_u \sum_j \alpha_u[j]\, N_{i,u,j}\;}_{\substack{\text{base load} \\ \text{(dwelling-count features)}}}
$$

where $N_{i,u,j}$ is the IPONDL dwelling count of level $j$ of variable $u$ in district $i$, and the $\alpha$'s are in MWh per dwelling.

The two blocks can use entirely different predictors — the heating correction keeps tenure / insulation / emitter; the base load takes household size, age, income, area class (household size is not even in the current candidate pool, and it is the canonical appliance-load driver). The group elastic net selects within each block naturally, and the whole thing remains one convex solve.

A side benefit for interpretation: `age_class` is currently the strongest electric-channel predictor, and it plausibly proxies base-load variation as much as heating behaviour — the decomposition would disentangle those two readings.

### Two honest caveats

1. **A new collinearity axis.** $E_i$ and $n_i$ are strongly correlated across districts; separating "scale the EPC mass" from "per-dwelling constant" leans on cross-district variation in conventional intensity $E_i / n_i$ and in composition. The elastic net's grouping handles the estimation side, but the identification-bounds analysis must be run on the joint design, and a soft trade-off ridge between $g_0$ and $\alpha_0$ is to be expected. An elegant mitigation is available: the SDES national non-heating statistic currently hard-coded per dwelling becomes a **single linear equality constraint** anchoring the base-load national total,
   $$
   \sum_i \hat b_i \;=\; \text{SDES national total},
   $$
   letting its spatial *distribution* be estimated while its level is pinned by the external statistic. That is strictly more information-respecting than the current constant, and it preserves convexity.
2. **Validation opportunity that doubles as a check.** The estimated heating/base decomposition of electricity can be compared against Enedis's modelled *thermosensible* split — which the paper currently avoids relying on. Agreement would be a strong face-validity result; disagreement would be genuinely informative either way.

The occupancy gate should apply to both components (vacant dwellings have no appliance load either), which the fixed-gate approach handles uniformly.

## 3. Recommended path

Frame it as a model-class upgrade with an ablation, not a silent swap:

1. **Parity check.** Fit the additive single-component model and confirm CV-RMSE parity with the multiplicative fit (expected, given moderate deviations — and if parity holds, the exact-convexity argument alone justifies the switch).
2. **Base-load block.** Add the second component, where the real accuracy gain is most likely on the electric channel (MAPE 14.6% with half the mass currently modeled as a national constant).
3. **Identification.** Rerun the identification bounds on the joint design with the SDES anchor.

The narrative also gets simpler, not more complicated: "additive main-effects correction, exact group elastic net, jointly convex, with a structurally separated base load" removes the two most attackable methods paragraphs (BCD approximation, floor-scaling) while keeping every interpretability claim intact — and it makes the MOB extension essentially free.

The one selling point given up is the multiplicative model's compounding semantics (effects as universal percentage scalings), which matches how the performance-gap literature reports things. If a reviewer pushes on that, the answer is that the two parameterizations agree to first order over the fitted range — demonstrable empirically with the parity check from step 1.
