# Meeting with Anders (Supervisor)

## Summary

### PDE solver discussion
- Question: is upwind (artificial diffusion) needed when the equation already has diffusion?
- Try substituting `w = exp(-w)` with a known solution
- Nonlinear Hamiltonian causes difficulties: instability, not finding solutions

### Plan B — Sensitivity / UQ on the ODE system
- If PDE proves too difficult, pivot to uncertainty quantification on the ODE model
- Study which parameters the solution is most sensitive to
- Ask the bank which parameters they have and which are most uncertain
- HJB known to be sensitive (per Anders)
- Start in the ODE system for simplicity, then verify against PDE
- Keep computation time low: solve on coarse grids, assume behavior transfers to finer grids

### Thesis structure
1. Model
2. ODE approximation — how to solve it
3. Sensitivity analysis on ODE parameters
4. PDE solver — compare against ODE in low dimension
5. Draw connection to the full problem

The thesis should be reproducible — give enough ODE background for the reader to recreate the work.

### Deadlines
- Anders away April 18–24
- Submit draft (even if incomplete) before April 18

---

## Raw notes

upwind -> artificiel diffusion, skulle det behövas då ekvationen har diffusion redan?

byta ut w = exp(-w) till en given lösning

svårigheter att lösa Hamiltonian som är icke linjär:
  - instabilt
  - inte hittar lösningen

ha någon plan B:
  - uncertainty quantificatino på ODE system
  - vilka parameterar är mest osäkra? Fråga banken
  - ODE system behöver data för parametrar
    - inte säkert att hinna använda riktig data
  -> behöver studie på hur osäkra olika parameterar är i ODE modellen

Fråga på banken va dodm har parameterar
hitta på hur lösningen beror på olika parameterar
 -> kan göra det mer eller mindre systematiskt 
 -> läsa rapporten

uncertainty quantification på ett grovt nät och se hur det påverkar där. 
om varje körning tar mycket tid så blir det inte lätt
  - lösa med stort numeriskt fel (kan alltid göra mindre)
  - anta att det beter sig likartat även fast man använder ett finare nät
  - viktigt att beräkningarna inte tar mycket tid (jätte besfärligt)
  - tar ett exempel på en parameter och bara visar att det är känsligt
    - HJB känt att vara känsligt (enligt the GOAT)
  - Var är mest känsligt? Vilka parametrar ska man studera?
  - börja i ODE systemet

studien skulle kunna göras i ODE systemet för enkelhetsskull

om man tar att ODE är känsligt så man kan göra det där

kan göra ett verifierande studie hur det stämmer överens

viktigaste med PDE lösningen är att kunna jämföra någon lösning i ODE lösning (låg dimension)

liknar mer rapporten med SIR modellen


prata med kollegor vad som kan vara känsligt / testa med körningar


rapportens structure:
modellen, ODE, hur man löser, känslighet, stämmer den bra överens med PDE lösningen, drar koppling till riktiga lösningen. 

man vill att den som läser rapporten ska kunna återskapa det jag har gjort. Bra att ge den bakgrund (ODE)


the GOAT är borta från 18e till 24e april
  - lämna in det jag har innan den 18e (tom den är ofärdig)