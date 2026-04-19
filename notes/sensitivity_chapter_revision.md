# Sensitivity Chapter Revision Notes

Organized into logical groups A--H. Tackle one group at a time.
Items marked with their original number from the review list.

---

## Group A: Factual Errors in Text vs Data

These are cases where the text states something that contradicts the actual Sobol index data. Must be fixed before submission.

### [x] A1. Summary table (Table 4.1) omits dominant parameters for QoI 6 [item 33]

The table caption says "dominant parameters (those with S\_Ti > 0.05)" but for QoI 6 (Net revenue) only lambda\_scale is listed. The actual data shows:
- gamma: S\_Ti = 0.069 (above 0.05)
- sigma\_EUR: S\_Ti = 0.054 (above 0.05)

**Fix:** Add gamma and sigma\_EUR to the QoI 6 row in Table 4.1. Or change the threshold in the caption (e.g., >0.10) and adjust all rows to match.

### [x] A2. Inconsistent "Interactions" column in Table 4.1 [item 34]

The sum-of-S\_Ti values and their labels:
| QoI | Sum S\_Ti | Label |
|-----|----------|-------|
| 1 | 0.994 | None |
| 2 | 1.054 | Small |
| 3 | 1.075 | Moderate |
| 4 | 1.261 | Significant |
| 5 | 1.047 | Small |
| 6 | 1.050 | None |

QoI 6 (1.050) is labeled "None" while QoI 5 (1.047) and QoI 2 (1.054) are labeled "Small". These are essentially the same value.

**Fix:** Either label QoI 6 as "Small" for consistency, or define explicit thresholds (e.g., <1.02 = None, 1.02--1.10 = Small, etc.) and apply uniformly.

### [x] A3. QoI 5 text says "No pair-specific demand parameter appears" [item 29]

The data shows beta\_1^EU has S\_Ti = 0.025 for QoI 5, which IS a pair-specific demand parameter and appears in the Sobol bar chart (threshold 0.02). The text also says the model is "nearly additive (sum S\_Ti = 1.05)" but then for QoI 3 (sum S\_Ti = 1.08) says "consistent with the presence of interactions." See also A5 below.

**Fix:** Either raise the bar chart threshold to 0.03 so beta\_1^EU drops out (and note it in text), or acknowledge beta\_1^EU as a minor contributor. More precisely: "No pair-specific parameter contributes more than 3%."

### [x] A4. QoI 6 text lists wrong parameter [item 31]

Text (sec 4.4.3): "several demand parameters from multiple pairs at the 1--2% level (alpha\_2^GU, beta\_2^EU, **beta\_2^GU**)."

The actual data: beta\_2^GU has S\_Ti = 0.011 (1.1%), which is at the noise floor. Meanwhile alpha\_1^GU has S\_i = 0.012 (same level) and is omitted. More importantly, **alpha\_2^EU** at S\_Ti = 0.044 (4.4%) is the third-largest contributor after lambda\_scale and gamma, and is mentioned separately in the "secondary" list, so this is OK. But the parenthetical list is misleading because beta\_2^GU barely qualifies at 1%.

**Fix:** Remove beta\_2^GU from the list (below noise floor). The parenthetical should be: "(alpha\_2^GU, beta\_2^EU)" at 2% each. Or raise the cutoff to 2%.

### [x] A5. Inconsistent "additive" vs "interacting" characterization [items 21, 28]

- QoI 2: sum S\_Ti = 1.054 -- text says "nearly additive"
- QoI 3: sum S\_Ti = 1.075 -- text says "consistent with the presence of interactions"
- QoI 5: sum S\_Ti = 1.047 -- text says "nearly additive"

The gap between 1.054 and 1.075 is only 0.02. Calling one "additive" and the other "interacting" overstates the distinction.

**Fix:** Use consistent language. Suggestion: anything below ~1.10 is "nearly additive", 1.10--1.20 is "moderate interactions", above 1.20 is "significant interactions" (only QoI 4 at 1.26). State the threshold once in the methodology section and apply uniformly. This also fixes A2.

### [x] A6. Heatmap description says "dark" but should say "light/yellow" [item 14]

Text (sec 4.3): "The GBP/USD and EUR/GBP pair-specific parameters (rows 11--20) are dark across all QoIs."

With the YlOrRd colormap, low values are light yellow, not dark. High values are dark red.

**Fix:** Change "dark" to "light yellow" or "near-zero" or simply "blank." The figure is being remade anyway (see Group D).

### [x] A7. Full data audit required [items 31, 33]

Given that multiple text-vs-data errors have been found, every numerical claim in sections 4.3--4.5 should be systematically verified against the printed Sobol indices. The QoI 1 section quotes S\_Ti values (72%, 15%, 9%, 3%) while using phrasing that implies S\_i ("accounts for X% of variance"). This works because QoI 1 is nearly additive so S\_i ~ S\_Ti, but is technically sloppy. Either quote S\_i values (72%, 17%, 11%, 4%) or clarify that the percentages are total-effect indices.

**Fix:** Go through each QoI section and cross-check every stated index value against the data tables in the notebook output. A spreadsheet-style check: for each QoI, list every parameter mentioned in the text and verify the S\_i and S\_Ti values match.

---

## Group B: Unit and Notation Consistency

### [x] B1. Hedge rate units: M$/day in Chapter 4 vs M$/s in Chapter 3 [item 43]

The ODE chapter (section 3.5.2, Figure 3.X) plots hedge rates in **M$/s** (the plotting function in `src/plotting.py:118` divides by DAY\_SECONDS = 86,400). The sensitivity chapter defines QoI 4 as M$/day and reports all values in M$/day (mean = -2,896 M$/day, CI = [-9,287, 0]).

Both are technically correct (the model uses days as the time unit, so the raw output is M$/day; the ODE figure converts to M$/s for display). But using different units across chapters is confusing.

**Fix options:**
1. **Convert sensitivity chapter to M$/s** (divide all hedge rate values by 86,400). Values become: mean ~ -0.034 M$/s, CI ~ [-0.107, 0] M$/s. These are small numbers but consistent with ODE figures.
2. **Convert ODE figure to M$/day** and unify everywhere on M$/day.
3. **Use both** with explicit conversion: define QoI in M$/day but note the M$/s equivalent for cross-reference to Chapter 3.

Recommendation: option 1 or 2. Pick one unit and use it everywhere. M$/day may be more intuitive for daily PnL context; M$/s matches the ODE figures. Either works, but pick one.

The net revenue (QoI 6) is in $/day, which should also be checked for consistency. The model chapter says "per unit time" without specifying; clarify.

### [x] B2. Sobol index format: % vs decimal [item 38]

In sections 4.3--4.4, some QoIs report Sobol indices as percentages ("72%", "6%") and others as decimals ("S\_i = 0.35, S\_Ti = 0.36"). This is inconsistent within the chapter.

**Fix:** Pick one format and use it consistently. Recommendation: use decimal notation (0.35) when giving both S\_i and S\_Ti, and percent (35%) when summarizing a single dominant contribution in prose. Whatever the choice, state it once and be consistent.

### [x] B3. "mean" vs "nominal" in forward UQ [item 18]

The KDE plots label the dashed red line as "Nominal" (the value at the paper's parameter set). The text sometimes says "mean" when referring to the sample average. These are different quantities (the nominal is one specific evaluation; the mean is the average over all 20,000 samples). The text correctly uses both, but should always be clear which is which.

**Fix:** No change needed in the text (the text does distinguish), but verify that every instance of "mean" refers to the sample mean and every instance of "nominal" refers to the paper's parameter values. In the forward UQ discussion: "The nominal value of 0.12 bps falls below the mean" -- this is correct usage.

### [x] B4. Net revenue sentence overflows page [item 32]

Text (sec 4.4.3): "The forward UQ gives a mean of \$502k/day, CV = 48%, and 90% CI = [\$124k, \$893k]/day."

This line overflows the page margin, probably because of the inline math.

**Fix:** Break the stats onto their own line or use a displayed equation. E.g.: "The forward UQ gives mean = \$502k/day, CV = 48%, 90% CI = [\$124k, \$893k]/day." or display as a small table.

---

## Group C: Remove/Rework Two-Currency References

A recurring pattern: the chapter repeatedly compares 3-currency results to a hypothetical 2-currency analysis. This confuses readers because the 2-currency case was never presented and isn't a natural comparison point.

### Items affected: 3, 9, 12, 20, 36

Specific locations:
- [x] **Intro** (line 8): "A two-currency analysis would miss all of these effects." [item 3]
- [x] **QoI 3 definition** (line 100): "This quantity is structurally zero in a two-currency model where GBP does not exist." [item 9]
- [x] **End of QoI definitions** (line 126): "QoIs 1, 2, and 4 have direct analogues in a two-currency setting. QoIs 3 and 5 probe mechanisms that exist only with three or more currencies." [item 12]
- [x] **QoI 3 results** (line 237): "This QoI probes a mechanism that does not exist in a two-currency model" [item 20]
- [x] **Discussion 4.5.1** (line 352): "The sensitivity rankings for QoIs 1, 2, and 4 reproduce what a two-currency analysis (with 8 parameters) would give." [item 36]

**Fix for all:** Remove or rework these sentences. The key insight to preserve is that cross-pair QoIs (3, 5) depend on rho and sigma\_GBP, which are unique to the multi-currency setting. State this positively: "QoIs 3 and 5 reveal sensitivity to rho and sigma\_GBP, parameters that have no effect on own-pair quantities." Don't frame it as "missing" from a 2-currency model.

For the intro [item 3], replace the sentence with motivation for *why* 3 currencies specifically: "Three currencies is the smallest setting that exhibits these multi-currency mechanisms while remaining computationally tractable for the 220,000-evaluation Sobol study. Adding more currencies would increase the parameter count but not introduce qualitatively new mechanisms." [item 2]

For the end-of-QoI-definitions paragraph [item 12], replace with: "This mix of own-pair (QoIs 1, 2, 4) and cross-pair (QoIs 3, 5) quantities, together with the aggregate (QoI 6), lets us test whether cross-currency coupling introduces new sensitivity dimensions beyond the own-pair parameters."

For the discussion paragraph [item 36], replace the 2-currency comparison with a more informative closing. For example, discuss what the sensitivity structure implies for model deployment at larger currency counts, or what the separation between microstructure and risk management parameters means practically.

---

## Group D: Figure Quality and Redesign

All figures in this chapter need to reach the same publication quality as the ODE chapter figures.

### [x] D1. Heatmap: colorbar overlaps group labels [item 15]

The colorbar scale on the right side of the heatmap overlaps with the parameter group labels (EUR/USD, GBP/USD). The group labels are placed at `N_QOIS + 0.3` but the colorbar `pad=0.15` isn't enough.

**Fix options:**
- Increase colorbar `pad` (e.g., 0.25)
- Move group labels to the left side of the heatmap
- Remove group labels and use horizontal separator lines (already present) with a legend
- Widen the figure (`figsize=(10, 9)`) to give more room

Also consider: the heatmap annotates values > 0.05 but the parameter table and text use different thresholds. Align them.

### [x] D2. Sobol bar charts: unexplained group coloring [item 16]

The per-QoI Sobol bar charts color bars by parameter group (blue for Global, red for EUR/USD, etc.) via the `GROUP_COLORS` dictionary. This coloring was inherited from another figure and is unexplained in the figure or caption. It adds visual noise without conveying information (the parameter names already indicate the group).

**Fix:** Remove the group coloring. Use a single color for all bars, with lighter shade for S\_i and darker for S\_Ti (the current legend already describes this). Match the style of the ODE chapter figures: clean, minimal, grayscale or single-accent-color.

Additional issues with the bar charts:
- The figure suptitle repeats the QoI name that's already in the x-axis label
- The legend shows generic gray patches that don't match the actual bar colors
- Consider: would a single combined figure (6 panels) be cleaner than 6 separate figures?

### [x] D3. Scatter plots: density hard to read [items 41, 44]

**eta vs hedge rate (Fig 4.X):** Most points cluster at hedge rate = 0 (inside dead zone), making the downward trend to -10,000 hard to see. The plot extends to -30,000 which wastes space.

**Fix suggestions:**
- Use a 2D histogram or hexbin instead of raw scatter for better density visualization
- Or use alpha-blending with a colormap (density-colored scatter)
- Clip y-axis to [-12,000, 500] to focus on the interesting range
- Add marginal distributions (histograms on axes)

**rho vs cross-skew (Fig 4.Y):** The red dashed line shows the nominal rho value, but in the eta scatter plot the red line shows the 1/eta reference curve. This inconsistency is confusing [item 44].

**Fix suggestions:**
- Remove the nominal rho vertical line (it doesn't add much)
- Instead, add a linear regression line to highlight the approximately linear relationship mentioned in the text
- Use consistent visual language: if red = reference curve, use it for the linear fit; if red = nominal, use it in both figures
- Same density visualization improvements as above

### [x] D4. Forward UQ density y-axis scales [item 26]

The KDE plots have very different y-axis scales: QoI 4 density peaks around 0.00025, QoI 5 around 80,000. This is expected (density integrates to 1, so narrow distributions have high peaks and wide distributions have low peaks). This is not an error.

**Fix:** Not strictly necessary, but adding the y-axis label "Probability density" makes it clearer. Could also normalize to unit peak height if the shape matters more than the scale, but standard KDE is fine for a thesis.

---

## Group E: Writing Improvements (Prose and Phrasing)

### [x] E1. "key strategy outputs" [item 1]

Intro: "...to obtain the probability distributions of the key strategy outputs..."

**Fix:** Replace with "...to obtain the probability distributions of each quantity of interest (QoI)..." since QoI is the term used throughout.

### [x] E2. Motivation for 3 currencies [item 2]

Intro: "We work with three currencies (USD, EUR, GBP), the smallest setting that includes..."

**Fix:** Add: "...the smallest setting that includes cross-pair correlations, cross-inventory effects on quoting, and cross-hedging. Using the minimal multi-currency configuration keeps the parameter count manageable (20 parameters vs. potentially 50+ for five currencies) while capturing every qualitative mechanism in the model."

### [x] E3. "exogenous" [item 6]

End of 4.1.1: "...the exogenous trade rates k\_i (negligibly small)."

**Fix:** Replace "exogenous" with "external" or "background." E.g., "the background trade rates k\_i (negligibly small)."

### [x] E4. Parameter names in table [item 7]

"Logistic shift" and "Logistic slope" are technical but acceptable. The readers of this thesis will have read Chapter 2 where the logistic model is introduced. Could add a brief parenthetical reminder.

**Fix (optional):** Rename to "Demand curve intercept" and "Demand curve steepness" if you want to be more intuitive. Or keep as is -- this is minor.

### [x] E5. Unnecessary sentence about determinism [item 8]

"All six are deterministic given the ODE solution and introduce no sampling noise."

This is true and relevant (it distinguishes the ODE-based SA from a Monte Carlo simulation-based SA where QoI evaluations themselves are noisy). But it reads as throat-clearing.

**Fix:** Remove, or shorten to a parenthetical: "...we extract six scalar QoIs (all deterministic given the ODE solution)."

### [x] E6. Missing equation reference for markup formula [item 17]

Text: "The optimal markup is therefore dominated by the maximization of f^n(delta) * delta over the logistic demand curve."

**Fix:** Add the equation reference: "The optimal markup~\eqref{eq:delta-star} is therefore..."

### [x] E7. "before the dead zone clips it" [item 27]

QoI 5: "the raw hedging momentum p\_EURGBP before the dead zone clips it."

**Fix:** Rephrase: "the raw hedging momentum p\_EURGBP, which measures the value function's incentive to cross-hedge before execution costs are applied."

### [x] E8. QoI numbers inconsistent in 4.5.1 [item 35]

The paragraph in 4.5.1 names some QoIs by description only and others by number:
"The tier spread differential is a microstructure quantity... The cross-pair quantities (QoIs 3 and 5)... Net revenue is dominated..."

**Fix:** Use both name and number consistently throughout: "The tier spread differential (QoI 1) is... The inventory skew (QoI 2) and hedge rate (QoI 4) are... The cross-pair quantities (QoIs 3 and 5)... Net revenue (QoI 6)..."

### [x] E9. Section 4.5.4 feels abrupt [item 47]

The PDE bridge section currently reads as an afterthought: "All results in this chapter are based on the ODE approximation... In Chapter 5, we develop a PDE solver..."

**Fix suggestions:**
- Open with a forward-looking framing: "The sensitivity analysis above rests on the ODE approximation. A natural question is whether the sensitivity structure persists when the exact Hamiltonians are used."
- Mention what specifically might differ: "The quadratic Taylor expansion is least accurate at large inventories and large hedging momenta -- precisely the regime where the hedge rate (QoI 4) is most uncertain."
- End with a concrete promise: "Chapter 5 addresses this by comparing ODE and PDE solutions at the extremes of the parameter ranges identified here."

---

## Group F: Inventory Position and QoI Design

Several items question whether 10 M$ is the right inventory level for the hedging QoIs.

### [ ] F1. QoI 4: hedge rate at y=10M$ EUR [items 10, 22, 23, 24]

At y\_EUR = 10 M$, the nominal hedge rate is -1,568 M$/day. The 90% CI reaches 0, meaning many parameter combinations land inside the dead zone (no hedging at all). This creates a heavily skewed distribution with a point mass near zero, making the QoI hard to interpret (CV = 111%).

At a larger inventory (e.g., y\_EUR = 20 M$), the hedge rate would be further from the dead zone edge, producing a more informative distribution. The sanity check cell in the notebook actually computes at y=20M$ and gets -8,137 M$/day (nominal, n=2000), which is well outside the dead zone.

**Arguments for keeping 10 M$:**
- 10 M$ is a realistic single-trade size (it's in the model's size ladder)
- The dead zone behavior IS the story: it shows that hedging is fragile near the boundary
- Changing the inventory would require rerunning the full 220,000-evaluation Sobol study

**Arguments for changing to 20 M$:**
- The distributions would be more informative (less point-mass at zero)
- Cross-hedging (QoI 5) might become active, revealing more structure
- The hedge rate CI would be more meaningful

**Recommendation:** Keep 10 M$ for the main analysis but consider adding a brief discussion (or appendix result) showing how the hedge rate sensitivity changes at 20 M$, to address the dead-zone boundary effect. This would also address items 23 and 24 about the cross-pair dead zone blocking cross-hedging.

### [ ] F2. QoI 5: cross-hedge momentum [item 11]

The cross-hedge momentum is studied instead of the actual cross-hedge rate because the dead zone blocks cross-hedging at y=10M$. At larger inventory (20M$), the actual cross-hedge rate might become nonzero for some parameter combinations, making it a viable QoI.

**Recommendation:** Same as F1. The momentum is a valid QoI that reveals the sensitivity structure of the cross-hedging incentive. But note in the text that at larger inventories, the dead zone would be breached and actual cross-hedging would occur.

### [ ] F3. Cross-pair dead zone and inventory dependence [item 23]

Text: "the EURGBP dead zone blocks cross-hedging for nearly all parameter combinations... so the EURUSD hedge rate is determined as if it were the only hedging channel."

This is partly an artifact of the chosen inventory level. At larger EUR inventories, the cross-hedge momentum exceeds the dead zone and EURGBP hedging activates.

**Fix:** Add a sentence: "At the chosen inventory of 10 M\$, the cross-hedge momentum (QoI 5) is smaller than the EURGBP dead zone for nearly all parameter combinations. At larger inventories, the dead zone would be breached and the EURUSD hedge rate would also depend on cross-pair parameters through the activated EURGBP channel."

### [ ] F4. Stressed inventory scenario [item 25]

Would a larger inventory "stress test" all QoIs? Quoting QoIs (1--3) would change because the momentum p increases with y. Hedging QoIs would change substantially. Revenue would change. This could be interesting follow-up work.

**Recommendation:** Not necessary for the current chapter, but worth mentioning in the discussion as a limitation/future direction: "The QoIs are evaluated at relatively modest inventory perturbations (10 M\$). At larger inventories, the quoting momentum deviates further from zero, potentially activating nonlinear effects in the demand curve that the sensitivity analysis at 10 M\$ does not capture."

---

## Group G: Structure and Repetition

### [x] G1. Lambda scale introduction order [item 4]

In 4.1.1, the lambda scale parameter appears in the parameter table (paragraph 1) but is only explained in paragraph 3 ("The arrival rates... are scaled by a single shared multiplier lambda\_scale...").

**Fix:** This is acceptable -- the table introduces the parameter and the explanation follows. But could move the explanation to paragraph 1 or add a brief note in the table caption: "The arrival rate scale lambda\_scale is a shared multiplier applied to all arrival rates (see text)."

### [x] G2. Repetition between intro and 4.1.1 [item 5]

The intro paragraph 1 explains which parameters are uncertain and why. Section 4.1.1 paragraph 2 repeats much of this.

**Fix:** Keep the detailed explanation in 4.1.1 and shorten the intro. The intro should say *that* the parameters are uncertain (one sentence); section 4.1.1 explains *why* and *how much*. Remove from the intro: "Exchange rate volatilities fluctuate... Client arrival rates depend... The logistic demand parameters... are fitted from historical fill-or-reject data..." -- move all of this to 4.1.1 if not already there.

### [x] G3. CV context [item 40]

The methodology section introduces CV but doesn't give ranges for what constitutes "low" or "high." Then the results section says "relatively robust" for CVs of 38--46%.

**Fix:** Add context in the methodology section, e.g.: "For reference, a CV below 25% indicates that the output is tightly constrained by the input uncertainty, 25--50% indicates moderate sensitivity, and above 100% signals that the output is dominated by parameter uncertainty." (Adjust thresholds as appropriate.) This gives the reader a frame before encountering the results.

### [x] G4. Table format inconsistency [item 39]

The summary table (4.5.1) lists QoIs without units; the forward UQ table (4.5.3) includes units. This is fine since they serve different purposes, but could be made more parallel.

**Fix (optional):** Add units to the summary table's QoI column. Minor improvement.

### [x] G5. Parameter reader burden [item 30]

20 parameters with Greek symbols across 3 pairs. The reader cannot remember them all. The parameter table helps, and the grouping (global vs pair-specific) helps.

**Fix (optional):** Consider adding a symbol glossary in the appendix, or a margin note at first use of each symbol in the results sections. This is standard thesis practice. Not urgent.

### [x] G6. ODE vs PDE model [item 48]

The sensitivity analysis is on the ODE model, and the text explains parameters in terms of ODE structures (matrix A, Riccati equation). The PDE model has the same parameters but different mechanics. The chapter doesn't discuss how parameters enter the PDE.

**Fix:** This is fine for the thesis structure. The sensitivity analysis IS on the ODE model. The PDE comparison in Chapter 5 validates whether the ODE sensitivity conclusions carry over. No change needed, but the bridge section (4.5.4) should be clearer about this (see E9).

### [x] G7. Net revenue CI reference [item 45]

Text (4.5.3): "The net revenue remains positive across the entire parameter space explored (90% CI = [$124k, $893k]/day)."

This value appears in the forward UQ table above but the text doesn't reference the table.

**Fix:** Add a reference: "...as confirmed by the forward UQ in Table~\ref{tab:sa-forward-uq}."

### [x] G8. "factor-of-two" hedge rate claim [item 42]

Text: "Underestimating market impact by a factor of two can lead to a factor-of-two overestimate of the optimal hedging rate."

Analysis: The hedge rate formula is xi* = (q - psi)/(2*eta). If eta halves (underestimate by 2x), xi* doubles (holding q constant). So the claim is correct for the direct effect of eta. However, the scatter plot shows substantial spread at each eta value (due to other parameters), so the *total* effect of a 2x eta error depends on the other parameters too. The statement is correct as a worst-case / direct-effect claim.

**Fix:** Slightly soften: "Underestimating market impact by a factor of two approximately doubles the prescribed hedging rate, as the formula xi* = (q - psi)/(2*eta) is inversely proportional to eta."

### [x] G9. More discussion like end of 4.5.3 [item 46]

The last paragraph of 4.5.3 (about the desk being profitable but not knowing how aggressively to hedge) is insightful. More discussion in this style would improve the chapter.

**Fix:** Add similar practical insights after each QoI group:
- After Set A: "The quoting policy is well-determined: a bank can set client spreads with moderate confidence even under parameter uncertainty, as long as the demand curves are estimated."
- After Set B: "The hedging prescription is the Achilles' heel: the desk knows it should hedge, and knows it will be profitable, but the optimal aggressiveness is highly uncertain."
- In the discussion: "For a bank deploying this model, the practical message is: invest calibration effort in eta and the demand curves, accept that hedging rates will need real-time adjustment, and trust that the quoting policy is robust."

---

## Group H: Notebook Fixes

### [x] H1. Missing `all_params` variable on computation path [item 49]

When running the computation cells (skipping "Load saved results"), the variable `all_params` is never defined. It is only created in the load cell as `all_params = np.concatenate([A_samples, B_samples], axis=0)`. The scatter plot cells use `all_params` and will crash with a NameError.

**Fix:** Add `all_params = np.concatenate([A, B], axis=0)` to the computation/save cell (after `fwd_qois = np.concatenate([f_A, f_B], axis=0)`).

### [x] H2. Sanity check cell uses wrong inventory level (BUG FOUND)

The sanity check cell (comparing n\_steps=500 vs n\_steps=2000) uses `y_long_eur = np.array([0.0, 20.0, 0.0])` for the high-accuracy reference, but `evaluate_qois()` uses `y_long_eur = np.array([0.0, 10.0, 0.0])`. This means the comparison is invalid: the 50--80% "relative errors" shown are NOT discretization errors but come from evaluating at different inventory levels.

Evidence: the tier spread differential (which doesn't depend on y\_long\_eur) shows zero error (2.64e-16), while all inventory-dependent QoIs show ~50% error, which is exactly what you'd expect from doubling the inventory (the skew scales roughly linearly with y).

**Fix:** Change the sanity check cell to use `y_long_eur = np.array([0.0, 10.0, 0.0])` and `y_long_gbp = np.array([0.0, 0.0, 10.0])` to match `evaluate_qois()`. The actual discretization error (n=500 vs n=2000) should be negligible.

---

## Summary: Priority Order

| Priority | Group | Effort | Impact |
|----------|-------|--------|--------|
| 1 | **A: Factual errors** | Medium | Critical -- wrong numbers in text |
| 2 | **B: Unit consistency** | Medium | Critical -- M$/day vs M$/s |
| 3 | **D: Figure quality** | High | High -- figures not publication-ready |
| 4 | **H: Notebook bugs** | Low | Medium -- reproducibility |
| 5 | **C: 2-ccy references** | Low | Medium -- reader confusion |
| 6 | **E: Prose fixes** | Low | Medium -- writing quality |
| 7 | **G: Structure** | Low-Med | Low-Med -- polish |
| 8 | **F: Inventory discussion** | Low | Low -- scope decision |
