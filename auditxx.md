# Populate generic estimate/lower_ci/upper_ci columns using reference q

Ha pensado durante 29s

I audited the current state of PR \#66, including the PR
description/commit history, the 126-file diff, the core statistical
code, divergence/bootstrap implementation, S4 integration, tests, and CI
configuration. The PR has grown to 170 commits and is merging devel into
stable; GitHub currently shows only one PR check, a CircleCI
build_and_upload check. G GitHub +1

Verdict

I would not merge this PR as-is.

The direction is substantially better than the earlier
implementation—especially removing ARIMA differencing from the
confirmatory q-functional analysis, fixing the divergence definition,
adding explicit pairing validation, and adding Monte Carlo tests. The PR
itself documents several important limitations. G GitHub +2

However, I found several issues that are serious enough to block a
statistical package release:

Severity Finding 🔴 P0 paired=FALSE can silently turn into
pair-respecting bootstrap merely because metadata contains a
pairing-like column 🔴 P0 Generic divergence estimate/lower_ci/upper_ci
columns silently use the nearest q to 1, even when q=1 was not requested
🔴 P0/P1 The PR contains multiple competing Tsallis-divergence
implementations with materially different semantics, making it possible
for R/C++/bootstrap paths to disagree 🟠 P1 The primary paired GAMM path
ignores supplied weights, while the unpaired path uses them; this makes
method selection change the estimand/variance model 🟠 P1 Pair
autodetection accepts arbitrary columns without validating that they
actually represent pairs 🟠 P1 Statistical validation is deliberately
skipped on Bioconductor, while the PR is intended to harden the
Bioconductor-facing package 🟠 P1 The PR is far too large/mixed for
reliable review and regression isolation

Below are the important details.

1.  🔴 paired=FALSE is not actually authoritative

This is probably my biggest API/semantic concern.

In .prepare_divergence_execution(), when bootstrap is enabled, the code
automatically searches for pairing metadata. If it finds something like
subject_id, it assigns pair_ids regardless of whether the caller
explicitly requested paired = FALSE. The only consequence of
paired=FALSE is a message:

“Using pair-respecting bootstrap resampling.”

The resulting pair_ids are still passed down into the computation. G
GitHub

The relevant logic is:

bootstrap enabled → detect pair IDs pair IDs found → use them
paired=FALSE only changes the informational message

That means:

calculate_divergence( se, bootstrap = TRUE, paired = FALSE )

can produce paired bootstrap confidence intervals if the object happens
to contain a subject_id, pair_id, etc.

That’s a dangerous API contract.

The documented parameter says paired controls the paired design, and the
implementation instead treats metadata autodetection as higher priority.
G GitHub

Why this matters

A subject_id column does not necessarily mean that the analysis is
paired.

For example, a dataset may contain:

subject_id = patient ID condition = disease/control

while the analysis intentionally treats groups independently because
subjects differ between groups, or because the particular analysis is
not paired.

The current implementation silently changes the bootstrap resampling
unit.

Recommendation

Make the precedence explicit:

if (paired) { detect / validate pair structure } else { pair_ids \<-
NULL }

If automatic pairing is desired, introduce a separate explicit option
such as:

paired = “auto”

rather than making FALSE mean “false unless metadata convinces us
otherwise.”

Blocker.

2.  🔴 Generic divergence estimates silently alias to another q

This is a concrete correctness bug.

.finalize_divergence_matrices() says:

# value (q=1)

q_ref \<- 1 q_idx \<- which.min(abs(q - q_ref))

and then uses that q to populate the generic columns. G GitHub

The problem is that which.min(abs(q - 1)) does not mean “q=1”.

For:

q = c(0, 0.5)

the generic estimate becomes the q=0.5 estimate.

For:

q = c(0, 2)

it becomes q=0 because of the tie.

So the returned object can effectively say:

estimate = …

while the user may reasonably interpret that as the canonical q=1
divergence.

This is especially problematic because the code explicitly describes q=1
as the reference q. G GitHub

Recommendation

Require exact q=1:

q_idx \<- which(abs(q - 1) \< tolerance)

if (length(q_idx) == 1) { … } else { estimate \<- NA_real\_ lower_ci \<-
NA_real\_ upper_ci \<- NA_real\_ }

Or remove the generic aliases entirely for multi-q results.

The S4 layer already demonstrates the safer philosophy: it resolves q
explicitly and errors rather than silently analyzing a different q. G
GitHub

Blocker.

3.  🔴 There are multiple Tsallis divergence implementations

This PR claims to have corrected the divergence implementation and says:

D_q(P\|\|Q) = (sum P_i^q Q_i^(1-q) - 1)/(q - 1) with no abs(). G GitHub

The main divergence implementation follows that formulation, including
explicit q=0 handling. G GitHub +1

But another implementation in bootstrap_divergence_helpers.R still
contains:

divergence \<- (1 - sum_term)/(q - 1) … divergence \<- abs(divergence)

G GitHub

So the repository contains two conceptually different implementations:

Main/C++ path (sum_term - 1) / (q - 1)

with no abs(). G GitHub

Helper path (1 - sum_term) / (q - 1) abs(…)

G GitHub

Even if these happen to yield equivalent positive values for the
intended domain in some cases, having separate mathematical definitions
in production code is unacceptable for a statistical package.

It creates exactly the sort of R-vs-C++ / point-estimate-vs-bootstrap
discrepancy that the PR is supposed to eliminate.

Recommendation

Create one authoritative scalar implementation and make every route use
it:

R scalar ↓ R vector ↓ C++ accelerated implementation ↓ bootstrap

Then add a test asserting equivalence across all paths for a grid such
as:

q = c(0, .1, .25, .5, .9, 1, 1.1, 2, 5)

with:

zeros disjoint support identical distributions highly skewed
distributions pseudocount 0 pseudocount \> 0

The existing mathematical-property test battery is a good foundation,
but the cross-implementation equivalence test is the missing piece. The
PR currently advertises a large divergence test suite, including
identity, non-negativity, KL limit, asymmetry, q=0 and support
properties. G GitHub +1

Blocker until the implementation paths are unified or explicitly proven
equivalent.

4.  🟠 Paired GAMM ignores supplied weights

The new paired nlme::lme path explicitly says:

“Weights are not used in this path.” G GitHub

Meanwhile the preprocessing code can generate heteroscedasticity-derived
weights or accept user-supplied weights. Those weights are passed to the
GAM/GAMM machinery for other paths. G GitHub +1

Thus model dispatch can change what the analysis does:

paired → nlme::lme → weights ignored unpaired → mgcv::gam → weights used

That isn’t merely an implementation detail. It can change the estimated
variance and consequently p-values.

The comments acknowledge the issue, saying the weights are “not directly
transferable to nlme varFunc classes.” G GitHub

Recommendation

Either:

implement the equivalent variance model in nlme, e.g. a suitable weights
= varFunc(…), or explicitly disable adaptive weighting for the
confirmatory paired path, with the API and result metadata saying so.

I would prefer option 2 unless there is a strong scientific reason for
weighting: the PR is already correctly recognizing that data-adaptive
weighting shouldn’t silently influence the confirmatory test. G GitHub

5.  🟠 Pair autodetection is too permissive

.detect_pair_ids() accepts the first existing column among:

paired_samples pair_id pair_samples subject_id patient_id

and returns it as a pairing structure as long as it has no NAs. G GitHub

But it doesn’t actually enforce the documented requirement that there be
at least two samples per pair.

So this is accepted:

subject_id A B C D E

with every subject appearing once.

That’s not a paired design.

The downstream validation only catches certain repeated-condition
structures; it doesn’t establish that a valid pair actually exists. G
GitHub

Recommendation

Validate:

at least two observations per pair; ideally exactly one control + one
treatment for complete pairs; explicitly classify incomplete pairs;
verify that the pair identifier is actually shared across the two
conditions.

And return NULL if none of the candidate columns passes those checks.

6.  🟠 The C++ “paired” bootstrap API itself assumes positional pairing

There is another architectural smell in the lower-level C++
implementation.

divergence_bootstrap_compute_cpp() says:

// Paired bootstrap requires equal-length vectors

and then treats element i of x and element i of y as a pair. G GitHub

That is a positional pairing assumption.

The newer flexible implementation supports explicit x_pair_ids and
y_pair_ids, which is much safer. G GitHub

So there are now multiple pairing APIs with different semantics:

paired=TRUE ↓ old C++ function: positional pairs

pair_ids supplied ↓ new flexible C++ function: explicit pair IDs

That is an API trap for future maintainers.

Recommendation

Deprecate/remove the positional paired implementation or make it
explicitly named something like:

paired_bootstrap_positional_cpp()

and make all production callers go through the explicit pair-ID
implementation.

7.  🟠 Bioconductor validation is being skipped

The PR advertises 22 Monte Carlo validation tests, but the NEWS/PR
description says these are skipped on Bioconductor builds via
skip_on_bioc(). G GitHub

Meanwhile the CI configuration adds a Bioconductor check step. G GitHub

This leaves a gap:

local / ordinary test environment ↓ 22 statistical validation tests

Bioconductor environment ↓ validation tests skipped

For ordinary package correctness this might be reasonable if the
simulations are too expensive. For a PR whose central claim is
statistical inference hardening, however, this weakens the strongest
evidence supporting the changes.

Recommendation

Split the suite into:

fast deterministic statistical invariants → must run everywhere
expensive Monte Carlo calibration → nightly/explicit validation job a
smaller representative Monte Carlo smoke test → run in PR CI

At minimum, the PR should publish the actual MC results used to justify
the claimed calibration.

8.  🟠 The paired correlation model is better, but the fallback
    semantics need stronger guarding

The new primary model uses:

corCAR1(~q \| subject/condition)

which is a good improvement for irregular q spacing. The code explicitly
explains why grid-index AR(1) is problematic for irregular q and gives
an MC example where the old structure has inflated type-I error. G
GitHub

The fallback is grid-index corAR1, with a warning if the q grid is
irregular. That’s good.

However, the fallback is still capable of producing a statistically
different model after a convergence failure.

That means:

same data same method different numerical behavior ↓ different
covariance model ↓ different p-value

This should be surfaced prominently in the result, not merely as
internal metadata.

The PR does record model_used, fallback_level, and
correlation_structure, which is a strong design decision. G GitHub

Recommendation

Make fallback visible in the primary result table and provide a strict
mode:

strict = TRUE

where failure of the prespecified covariance structure returns NA/error
instead of silently changing the inferential model.

9.  🟡 The GAM slope effect size is not directly comparable with GEE

The code itself acknowledges this:

GAM computes difference in predicted slopes; GEE extracts a q×group
coefficient. G GitHub

This is correct to document, but the package needs to be very careful
not to expose these under a common slope_diff column and encourage
downstream comparisons.

For GAM:

ΔH/Δq over \[min(q), max(q)\]

For GEE:

coefficient of q × condition

Those are only equivalent under a linear specification.

Recommendation

Rename them:

slope_diff_predicted interaction_coefficient

or at least include an effect_size_definition field.

10. 🟡 The q-grid spline choice is defensible but deserves a stronger
    sensitivity test

The new spline basis uses:

k = max(min_k, min(max_k, n_q_unique - 1))

and explicitly avoids trying to infer curve complexity from CV. G GitHub
+1

That’s much cleaner than the previous adaptive approach.

But it is still a modeling choice that directly affects the interaction
test.

The PR has q-grid invariance tests, which is excellent. G GitHub

I’d add:

df = 2, 3, 4, 5

sensitivity tests on the same simulated curves, especially when there
are only 4–6 q values.

11. 🟡 The PR’s scope is itself a major risk

This is not a style complaint.

The PR currently represents:

170 commits 126 changed files version changes statistical-model changes
divergence rewrite C++ changes S4 changes plotting changes performance
refactoring documentation CI changes hundreds of tests multiple merges
from stable into devel. G GitHub +1

The commit history includes many “Fix tests”, “Fix issues”, “Fix bugs”,
“Increase coverage”, “Merge branch stable”, etc. G GitHub

This makes it exceptionally hard to establish that:

“the code in this PR corresponds exactly to the statistical claims in
the PR description.”

I’d strongly recommend splitting this into at least:

Divergence mathematical correction SAIT statistical-model changes S4/API
hardening Performance refactor CI/documentation

The statistical changes should be independently reviewable.

What I think is good

There is quite a lot here that is genuinely strong.

Removing ARIMA differencing

This is conceptually much cleaner. Treating q as a deterministic
functional argument rather than pretending it is a temporal index avoids
testing a derivative of H(q) when the scientific question is about the
original functional interaction. The PR explicitly documents this
rationale. G GitHub

Continuous-q correlation

Moving from grid-position AR(1) toward corCAR1(~q \| subject/condition)
is a very good correction for irregular q grids. G GitHub

Pair validation

The new validation rejecting ambiguous repeated observations within a
pair is good defensive programming. G GitHub

Explicit post-selection separation

The LASSO/ElasticNet result being explicitly exploratory and prohibited
from modifying the confirmatory model is exactly the right conceptual
boundary. G GitHub

S4 staleness metadata

Flagging inferential results after subsetting rather than pretending
they remain valid is particularly good. G GitHub

Plotting semantics

The change away from presenting averaged gene-wise CI bounds as a global
bootstrap CI is a meaningful statistical improvement. G GitHub

Merge recommendation

Request changes.

My minimum merge gate would be:

Must fix paired=FALSE must never activate pair bootstrap. Never alias
generic estimate to the nearest q. Unify the Tsallis divergence
implementation across R/C++/bootstrap. Add cross-path equivalence tests
for all q edge cases. Make fallback covariance/model changes explicit in
the user-visible result. Strongly recommended Validate pair columns
structurally before accepting them. Decide explicitly what happens to
heteroscedastic weights in paired models. Run at least a representative
statistical calibration suite in CI. Split the PR into independently
reviewable statistical/API/performance changes. Overall score

I’d rate the technical direction ~8/10, but the release readiness ~4/10.

The reason for the gap is that this is a statistical package: a handful
of silent semantic mismatches are considerably more dangerous than
ordinary software bugs because the code can produce perfectly
plausible-looking p-values and confidence intervals while answering a
different question.

The two issues I’d treat as absolute blockers are the paired=FALSE
override and the silent nearest-q aliasing. The multiple divergence
implementations are the next major concern because they undermine
confidence that point estimates and uncertainty estimates are actually
based on one mathematical definition. G GitHub +1

G G Fuentes
