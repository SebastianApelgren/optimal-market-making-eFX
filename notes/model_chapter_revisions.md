# Model Chapter (ch:model) — Revision Notes

Feedback from review, organized by priority and section. Each item has a status marker:
- [ ] = not started
- [x] = done

---

## Group A: Notation issues (1, 3, 6)

All three concern the `|_mid`, `|_hedge`, `|_cont` subscript notation on `dY`.

**Problem:** The subscripts `dY|_mid`, `dY|_hedge`, `dY|_cont` are informal shorthand. The original paper does not use them. They look like conditioning notation to some readers.

**Fix:** Replace with explicit prose. Instead of `dY_i|_mid = ...`, write the inventory SDE in full as a single equation that combines all three sources (mid-price, hedging, jumps), then refer to individual terms by name in the text. Alternatively, use the phrasing "the contribution from mid-price movements is..." without the subscript notation.

- [ ] Rewrite eq:inventory-mid (Section 2.1) — drop `|_mid`, present as "between trades, the inventory evolves as..."
- [ ] Rewrite the hedging drift display (Section 2.3) — drop `|_hedge`, say "hedging contributes an additional drift"
- [ ] Rewrite the Ito expansion (Section 2.5) — drop `|_cont`, present as "the continuous component of the expected change in theta"

---

## Group B: Factual checks (2, 4)

### Item 2: Arrival rates same across tiers?
**Location:** Section 2.2, last paragraph.
**Claim:** "the arrival rates are the same for both directions of a given pair and for all tiers."
**Status:** Correct per the code — `lambdas_per_day` is one array per canonical pair, shared by all tiers. Each tier contributes independently to the Hamiltonian with the same lambda. The tiers differ only in (alpha, beta).
**Action:** The sentence is accurate but could be clearer. Rephrase to: "the arrival rates depend only on the trade size and the pair, not on the tier or direction. The tiers differ solely in their demand-curve parameters."

- [ ] Rephrase arrival rate sentence for clarity

### Item 4: k_i ~ 10^-7 claim
**Location:** Section 2.3, last paragraph.
**Status:** Correct. Code: `k_EUR = 5e-3 * BP = 5e-7`. The code does NOT set k=0; it computes the impact terms but they are negligibly small. Statement is accurate.
**Action:** No change needed, but could add "and are retained in the numerical implementation" if desired.

- [ ] Optionally clarify that k_i are used in code but negligible

---

## Group C: Phrasing and word choice (5, 7, 8, 10, 11, 13)

### Item 5: Sum over fills notation
**Location:** Section 2.4, eq:objective.
**Problem:** `\sum_{\text{fills in}(t,T]}` is informal. The clarifying sentence below helps, but the notation itself is unusual.
**Options:**
(a) Keep as is — the clarifying sentence makes it clear
(b) Replace with a counting-process integral: `\int_t^T \sum_{(i,j),n,k} z_k \delta dN^n_{ij,k}` where N is the counting process
(c) Write as expected rate: `\int_t^T \sum ... \lambda f(\delta) z_k \delta \, ds` (but this anticipates the HJB)
**Recommendation:** Option (a) is fine for readability; option (b) is more formal but heavier. User to decide.

- [ ] Decide on notation for the sum over fills

### Item 7: "superimposed"
**Location:** Section 2.5, Client trades paragraph.
**Fix:** Replace "superimposed on" with "in addition to" or "alongside." E.g., "Client trades arrive as Poisson jumps alongside the continuous dynamics."

- [ ] Replace "superimposed"

### Item 8: "summand"
**Location:** Section 2.5, two instances.
**"summand" = "a term in a sum"** — standard math vocabulary, but if it feels unnatural:
- "Each term in the sum combines..." instead of "Each summand combines..."
- "...enters only through its own term in the sum, so..." instead of "...enters only its own summand, so..."

- [ ] Replace "summand" (two instances)

### Item 10: Bracket inconsistency
**Location:** Section 2.5, the rearrangement equation.
**Problem:** Left side uses `[...]`, right side uses `(...)`. Looks like the round brackets might be function arguments.
**Fix:** Use square brackets on both sides, or round brackets on both sides.

- [ ] Make brackets consistent in the rearrangement equation

### Item 11: Convex duality sentence
**Location:** Section 2.5, after eq:hedging-hamiltonian.
**Current:** "a connection to convex duality that will yield the closed-form evaluation in Section 2.6."
**Simpler:** "which we evaluate in closed form in Section 2.6." (Cut the convex duality reference — it's a nice observation but not essential here.)

- [ ] Simplify the Legendre-Fenchel sentence

### Item 13: "infinitesimal generator" phrasing
**Location:** Section 2.5, discussion of the four HJB terms.
**Current:** "it is the infinitesimal generator of the log-normal inventory process"
**Simpler:** "it arises from applying Ito's formula to the inventory process (2.2)" or "it reflects how the inventory's random fluctuations affect the value function."

- [ ] Simplify the diffusion term description

---

## Group D: Content/structural questions (9, 12, 14, 15, 16, 17)

### Item 9: DPP explanation too terse
**Location:** Section 2.5, "Assembling the HJB equation" paragraph.
**Problem:** The sentence "the expected infinitesimal change in theta plus the flow payoff equals zero" is hard to parse without DPP background.
**Options:**
(a) Add 2-3 sentences of DPP intuition before the derivation. E.g.: "If the dealer is already following the optimal strategy, then no short-term deviation should be able to improve the expected payoff. Over [t, t+dt], the dealer collects flow income (client revenue minus hedging costs minus risk penalty) and the value function changes due to the inventory dynamics. Optimality requires that these two effects balance exactly — if the sum were positive, the strategy could be improved by delaying, and if negative, by acting earlier."
(b) Keep it brief but rephrase: "Optimality requires that the value function cannot be locally improved: the expected change in theta from the inventory dynamics, plus the instantaneous payoff, must equal zero."
(c) Add a brief aside or footnote explaining the DPP.
**Recommendation:** Option (a) adds the most clarity with minimal length. User to choose.

- [ ] Expand DPP motivation (choose option a, b, or c)

### Item 12: Q and H notation
**Location:** Section 2.5, eq:hjb.
**Status:** The PDE chapter already uses Q(y,theta) and H_hedge(y, nabla theta). The paper uses similar notation. Keeping this is consistent.
**Action:** No change needed. Q for "quoting" and calligraphic H for "hedging Hamiltonian" is clear.

- [ ] No change (confirmed consistent with PDE chapter)

### Item 14: Curse of dimensionality — how does complexity grow?
**Location:** End of Section 2.5.
**Current:** Just says "infeasible due to the curse of dimensionality."
**Fix:** Add a sentence like: "The number of grid points grows exponentially with d: a grid with n points per axis requires n^d points in total. For the resolution used in Chapter 4 (n = 301), the two-currency case has roughly 9 x 10^4 points, but a five-currency case would require 301^5 ~ 2.5 x 10^12 — far beyond what fits in memory." This makes the exponential growth concrete.

- [ ] Add concrete dimensionality scaling example

### Item 15: Are closed-form Hamiltonians used in the ODE chapter?
**Status:** The ODE approximation uses the Taylor coefficients alpha_0, alpha_1, alpha_2 of H(p) around p=0. These are computed numerically from H(p) (which in turn uses the Lambert W formula). So the ODE chapter needs H(p) but uses it indirectly — it evaluates H and H' and H'' at p=0 to get the Taylor coefficients. The Lambert W closed form is the mechanism, but the ODE chapter's key object is the quadratic approximation, not the exact Hamiltonian.
**Action:** Rephrase the intro to Section 2.6 to say "These closed forms are used directly by the PDE solver (Chapter 4) and indirectly by the ODE approximation (Chapter 3), which takes a Taylor expansion of the quoting Hamiltonian around zero momentum."

- [ ] Rephrase Section 2.6 intro re ODE vs PDE usage

### Items 16 & 17: Should Section 2.6 live in ch:model or elsewhere?
**The argument for keeping it in ch:model:**
- H^n(p) and H_ij(q) are defined in the HJB (eq:hjb). Their closed forms are properties of the model choices (logistic f, quadratic+linear L), not numerical methods.
- Both ch:ode and ch:pde reference these formulas. Putting them in ch:model avoids duplication and makes both later chapters lighter.
- The Lambert W derivation is short (~15 lines of math) and flows naturally from the Hamiltonian definition.

**The compromise (recommended):**
- Keep the mathematical content (FOC, Lambert W, hedging three-case analysis) in ch:model.
- Remove the forward references to numerical stiffness from Section 2.6.2. Specifically, cut the sentence about "dominant source of numerical stiffness when solving the HJB equation on a grid, as we discuss in Chapter 4." Replace with a purely mathematical observation: "For the paper's parameters, eta ~ 10^-9, so 1/(2eta) ~ 5 x 10^8. Even a modest value function gradient outside the dead zone produces a large optimal hedge rate."
- Similarly, cut "This observation foreshadows a key finding of the sensitivity analysis in Chapter 5" from 2.6.1. Just state the mathematical fact.

- [ ] Remove forward references to numerical stiffness from 2.6.2
- [ ] Remove sensitivity foreshadowing from 2.6.1
- [ ] Keep the Lambert W and hedging derivations in ch:model (mathematical content stays)
