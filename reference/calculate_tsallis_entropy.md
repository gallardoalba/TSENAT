# Calculate Tsallis entropy for a vector of transcript-level expression values of one gene.

Calculate Tsallis entropy for a vector of transcript-level expression
values of one gene.

## Usage

``` r
calculate_tsallis_entropy(
  x,
  q = 2,
  norm = TRUE,
  what = c("S", "D", "both"),
  log_base = exp(1)
)
```

## Arguments

- x:

  Vector of (non-negative) expression values.

- q:

  Tsallis entropy parameter (q \> 0). Can be a single value or a numeric
  vector. Default is 2.

- norm:

  If `TRUE`, the entropy values are normalized by their theoretical
  maximum for the number of transcripts (so values lie in \[0,1\]).

- what:

  Which quantity to return: `"S"` for Tsallis entropy (S_q), `"D"` for
  Hill numbers (D_q), or `"both"` for a list with both.

- log_base:

  Base of the logarithm used for Shannon limits and normalization
  (default: `exp(1)`).

## Value

A numeric vector (named when length(q) \> 1) for `what = "S"` or `"D"`,
or a list with components `$S` and `$D` when `what = "both"`.

## Details

Implements S_q = (1 - sum p^q) / (q - 1) and D_q = (sum p^q)^(1/(1-q))
with the q-\>1 limits given by Shannon entropy and the exponential of
Shannon respectively. Natural logarithms are used for the q-\>1 limit
and normalization.
