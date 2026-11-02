# Approximate Poisson-leap acceptance, before the first approximate run

Algorithm: adaptive explicit Poisson leaps, halve and redraw a rejected negative
proposal, with an exact Direct event after 20 rejections. This is explicitly
approximate and rejection can add bias. It is not advertised as the full CGP
critical-reaction algorithm. Align steps to save times and use checked atomic
population updates; never clamp negative populations.

Fixed birth/death control: initial 30, immigration 20, death 0.5, horizon 2;
10000 paths each at epsilon 0.1, 0.05, 0.025; independent stream roots 9101..9103.
The exact law is Binomial(30, exp(-1)) + Poisson(40*(1-exp(-1))).
Family alpha 0.001 over three CDF rows: DKW epsilon_stat=sqrt(log(6/alpha)/(2n)).
Allow approximation CDF discrepancy <= epsilon/2 + epsilon_stat. Mean discrepancy
<= epsilon*exact_mean + 4.5*sqrt(exact_variance/n). Report raw discrepancies,
all failures and non-negative/completion accounting; do not require monotonic
sample errors. These are declared envelopes for this fixture, not universal bounds.

Nonnegativity, population/event/leap/output limits, strict duplicate/unknown fields,
zero rates and cancelled work have deterministic regression checks. Incomplete
ensembles have no unconditional statistics. Approximate output has its own schema
and assessment; legacy exact result/archive contracts remain unchanged.
