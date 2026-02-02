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

  Tsallis entropy parameter (q \> 0). Scalar or numeric vector (default:
  2).

- norm:

  Logical; if TRUE, normalize entropy by its theoretical maximum (values
  in \[0,1\]).

- what:

  Which quantity to return: "S" (Tsallis entropy), "D" (Hill numbers),
  or "both".

- log_base:

  Base of the logarithm used for Shannon limits and normalization
  (default: `exp(1)`).

## Value

For \`what = "S"\` or \`what = "D"\`: a numeric vector (named when
length(q) \> 1). For \`what = "both"\`: a list with components \`\$S\`
and \`\$D\`.

## Details

Implements S_q = (1 - sum p^q)/(q - 1) and D_q = (sum p^q)^(1/(1-q)).
Uses the q-\>1 limits (Shannon entropy and its exponential). Natural
logarithms are used for the q-\>1 limit and for normalization.

## Examples

``` r
# Basic usage with a small numeric vector
x <- c(10, 5, 0)
calculate_tsallis_entropy(x, q = c(0.5, 1, 2), norm = TRUE)
#>     q=0.5       q=1       q=2 
#> 0.5380048 0.5793802 0.6666667 
```
