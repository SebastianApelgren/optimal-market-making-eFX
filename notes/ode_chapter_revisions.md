# ODE Chapter (Ch. 3) — Revision Items

Reviewed 2026-04-14. Items grouped so each group can be tackled in one focused session.

---

## Group A — Wording and clarity in sections 3.1–3.2 (intro through ansatz)

- [x] **1. "nonlocal coupling" in the intro (line 4)**
What does "nonlocal" mean in "the nonlocal coupling through the quoting Hamiltonian"?
It refers to the fact that the quoting Hamiltonian evaluates θ at *shifted* points θ(y + z·d_{ij}), not just at y — so it couples each grid point to its neighbours at distance z. This is explained in the model chapter (line 193 of model.tex) but never clarified in the ODE intro. Either add a brief parenthetical or remove the jargon.

- [x] **2. PDE equation reference in the intro (line 4)**
"...make the PDE challenging to solve numerically, as we shall see in Chapter 5" — should include a direct reference to the PDE equation (e.g. `\eqref{eq:hjb}`) so the reader can look at what is meant.

- [x] **3. Alpha_1 and alpha_2 descriptions (section 3.1, lines 30–32)**
The descriptions are somewhat hard to follow:
- For α₁: the envelope theorem sentence is dense. Could simplify to something like: "Increasing the inventory cost p by dp lowers the net payoff δ − p by dp, so the expected loss at the optimum is f*(0)·dp."
- For α₂: "It controls how quickly the optimal profit deteriorates as the momentum departs from zero" — could be phrased more intuitively, e.g. "It measures the sensitivity of the quoting profit to inventory pressure."
The paper is similarly brief, so this is about our own exposition, not about matching the paper.

- [x] **4. Setting ψ = 0 in the hedging Hamiltonian (end of section 3.1, line 34)**
Two parts that need clarification:
- "set ψ_{ij} = 0 in the hedging Hamiltonian": explain *why* — the proportional cost ψ creates a dead zone (kink at q = 0) that makes H non-differentiable, so a Taylor expansion around zero wouldn't work. Setting ψ = 0 makes H purely quadratic (= q²/4η).
- "The proportional cost is retained in the formula for the optimal hedge rate after the value function has been computed; it only drops out of the Riccati ODE derivation": this means ψ is ignored *only* when deriving the ODE for A and B, but once A and B are known, the actual hedge-rate formula uses the full piecewise rule with the dead zone |q| ≤ ψ. Could be worded more plainly.

- [x] **5. "the algebra" (section 3.2, line 53)**
"The detailed algebra of this coefficient matching" — consider replacing with "the detailed derivation" or "the coefficient matching calculation."

- [x] **6. Where does the θ̂ ansatz come from? (section 3.2, line 42–46)**
The equation θ̂ = −y^T A(t) y − y^T B(t) − C(t) is introduced as "we seek an approximate value function of the form" without motivation. Add a sentence explaining that this is the most general quadratic form in y (with time-dependent coefficients), and that a quadratic value function is the natural companion to a quadratic Hamiltonian approximation — substituting both into the HJB makes everything polynomial in y.

- [x] **7. "This coupling is what links hedging to the quoting strategy through the matrix A" (end of section 3.2, line 97)**
"Coupling" is vague. Simplify to something like: "This is how the hedging cost structure feeds back into the quoting strategy: the hedging terms amplify the A² coefficient in the Riccati equation, making A (and hence the markup sensitivity to inventory) larger."

---

## Group B — Wording and clarity in section 3.3 (policy formulas)

- [x] **8. Section 3.3 overall tone**
The section uses dense/formal language. Could be made more readable by, e.g., leading each subsection with a one-sentence plain-English summary before the formula. Review for specific word choices that feel overly formal.

- [x] **9. "the proportional cost ψ_{ij} reappears here" (section 3.3.2, line 126)**
The displayed equation for q_{ij} (eq. 3.15) does *not* contain ψ. The ψ appears in the *next* step: the piecewise optimal-hedge-rate formula ξ* (referenced as eq. from the model chapter), which has the dead zone |q| ≤ ψ. The text should say explicitly that ψ enters through the *hedge rate formula applied to q*, not through q itself.

---

## Group C — Section 3.4 (parameters): text, scope, and formatting

- [x] **10. "for the remainder of the thesis" scope claim (line 133)**
The text says "we use the parameter set from Table 1 of [5]" for the whole thesis, but:
- The PDE chapter uses only 2 currencies (computational constraint).
- The sensitivity analysis also uses 2 currency pairs.
So the *parameter values* (σ, α, β, η, etc.) are the same throughout, but the *number of currencies* varies. The sentence should be qualified: e.g. "we use the parameter values from Table 1 of [5]; later chapters restrict to two or three currencies for computational reasons while retaining the same calibrated values."
Discussion: Should the sensitivity analysis consider more currency pairs? Should we remind the reader why we use 5 currencies here? (Answer: to reproduce and validate against the paper's results; the 5-currency case is the paper's main example.)

- [x] **11. Table formatting (line 139–173)**
The wider table (Table 3.1 for direct pairs) overflows the right margin. Fix the LaTeX layout — options: `\small` or `\footnotesize` font, `tabularx`, or `\resizebox`.

- [x] **12. Execution cost parameter name (line 135)**
"The execution costs are smallest for the most liquid pair (EURUSD)" — specify which parameter: e.g. "The execution costs (ψ and η) are smallest..."

- [x] **13. T = 0.05 days motivation (line 137)**
The time horizon T = 0.05 days (72 minutes) is used without explanation. In the paper, this is chosen because the strategy converges to its stationary value well before t = 0 — so T is effectively "long enough for stationarity." The text on line 99 says this implicitly ("the solution converges to a stationary value well before t = 0") but the parameter section should briefly motivate the specific choice, e.g. "The horizon T = 0.05 days is chosen to be large enough that the Riccati solution has converged to its stationary limit."

- [x] **14. "dollar-based legs" (line 159)**
"The correlation coefficient ρ describes the correlation between the two dollar-based legs" — this means the correlation between the two exchange rates X/USD and Y/USD that make up the cross pair. Rephrase to e.g.: "...the correlation between the constituent USD-denominated exchange rates."

---

## Group D — Section 3.5.1: bid/ask direction (potential error)

This group is the most critical — it involves a possible error in either the text or the code.

- [x] **15. Is the bid-ask distance really the "spread"? (line 187)**
Yes, the vertical distance between bid and ask markups is the full spread. The markup δ is the amount the dealer adds on top of the inventory cost: the bid *quote* is mid − δ_bid and the ask *quote* is mid + δ_ask. The plot shows δ_bid (positive) and −δ_ask (negative), so the vertical gap = δ_bid + δ_ask = the spread. (This is just a clarification for the reader, not a fix needed.)

- [x] **15b. Figure plotting convention — change bid/ask sign convention**
The current code plots bid as +δ_bid (positive) and ask as −δ_ask (negative). This is confusing because the negation on the ask has no economic meaning — it's just a display choice. Change to plot the quote relative to mid:
- Ask = +δ_ask (positive, above mid)
- Bid = −δ_bid (negative, below mid)
This way the Y-axis label "Bid and Ask Quotes (bps)" means what it says, ask is naturally on top, bid on the bottom, and the vertical gap is the spread. Both `plotting.py` and `ode_chapter_figures.py` need updating.

- [x] **16. Bid/ask skew direction — TEXT IS WRONG (lines 189)**
The text says: "As the GBP inventory increases, the dealer *raises* the ask markup (making it more expensive for clients to buy GBP) and *lowers* the bid markup (making it cheaper for clients to sell GBP)."

This is backwards. The correct economics:
- Long GBP → want to sell GBP → **lower ask markup** (cheaper for clients to buy GBP from dealer → attracts outflow)
- Long GBP → don't want more GBP → **raise bid markup** (worse bid price for clients → discourages clients selling GBP to dealer)

Code analysis confirms the code is correct:
- `markup(ccy_pay=GBP, ccy_sell=USD)` = bid markup. Momentum p increases with y_GBP → δ* increases → bid markup rises with inventory. ✓
- `markup(ccy_pay=USD, ccy_sell=GBP)` = ask-side markup. Momentum p decreases with y_GBP → δ* decreases → ask markup falls. ✓

**The text description should be corrected to match the code and the paper's Figure 1.** The plot itself is correct; only the verbal description is wrong.

- [x] **17. Review rest of section 3.5.1 analysis**
Given the error in item 16, re-read the EURUSD and EURGBP paragraphs (lines 191–194) carefully and make sure the skew direction descriptions are consistent with the corrected convention. The economic logic in these paragraphs may also need the directions flipped.

---

## Group E — Section 3.5.3: inventory distribution analysis and figure

- [x] **18. "superimposed" + complex phrasing (line 229)**
"The contours of the risk function γ/2 y^T Σ y are superimposed: the inventory distribution is clearly shaped by these risk contours, with the principal axis of the distribution aligned with the direction of minimum risk."
Rephrase, e.g.: "The dashed contours show levels of the portfolio risk γ/2 y^T Σ y. The inventory distribution follows the shape of these contours: it extends further in the direction where the risk is smallest."

- [x] **19. Complex follow-up sentence (line 229–230)**
"This alignment arises because the dealer penalizes risk proportionally to y^T Σ y through the mean-variance term in the objective, and the optimal strategy acts to keep the inventory in the low-risk region."
Simplify, e.g.: "This is because the dealer's objective penalises portfolio risk, so the optimal strategy keeps inventories close to the low-risk region."

- [x] **20. Inventory distribution figure range (figure_inventory_distribution in ode_chapter_figures.py)**
The current range is set by `max(99.5th percentile, 50)` and ends up being very wide relative to where the density sits, making the distribution hard to see. Change the code so the range is tighter (e.g. ±20 or ±30) to zoom in on the density. Keep this as a notebook-runnable parameter so the final range can be eyeballed interactively.

- [x] **21. Correlation direction — TEXT IS WRONG (line 231)**
"The dealer tolerates larger combined EUR-GBP positions when they are positively correlated (both long or both short) than when they are negatively correlated (long one, short the other), because the portfolio variance is lower in the former case."

This is backwards. For EUR and GBP with *positive* correlation ρ > 0:
- y^T Σ y includes the cross-term 2ρ·σ_EUR·σ_GBP·y_EUR·y_GBP
- When both long (y_EUR > 0, y_GBP > 0): cross-term is **positive** → variance is **higher**
- When one long, one short (y_EUR · y_GBP < 0): cross-term is **negative** → variance is **lower**

So the dealer tolerates larger positions when they have **opposite signs** (one long, one short), because the positions partially hedge each other — the portfolio variance is lower. The text should be corrected to say: "The dealer tolerates larger positions when the EUR and GBP inventories have opposite signs (long one, short the other), because the positive correlation between the exchange rates means such positions partially offset each other, reducing portfolio variance."

The distribution should be elongated along the y_EUR ≈ −y_GBP direction (anti-diagonal), not along y_EUR ≈ +y_GBP. Verify this matches what the current figure actually shows.
