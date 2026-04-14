# Table Rendering Issue: Root Cause & Solution

## Problem

Tables in `articles/TSENAT.md` were not expanding horizontally in the
generated HTML (`TSENAT.html`). Despite CSS styling attributes
attempting to set `width: 100%`, the tables remained constrained to
narrow widths, making them difficult to read.

### Affected Tables

1.  **Table 3**: Pairwise Tsallis divergence estimates (4 rows × 4
    columns)
2.  **Table 4**: Top 6 genes by effect size (6 rows × 7 columns)

## Root Cause

The CSS styling attributes were **incorrectly placed on the caption
paragraph** instead of the **table element itself**.

### Original Markdown Structure

``` markdown
|         |  q_0.01 |  q_0.05 |   q_0.1 |
|:--------|--------:|--------:|--------:|
| FOXJ2   | 0.04377 | 0.04633 | 0.04952 |
...

**Table 3 \| Pairwise Tsallis divergence estimates.** Divergence
(distance) between conditions... {.table .table .table-striped
.table-hover style="width: 100%; margin-left: auto; margin-right: auto;"}
```

**Problem**: The Pandoc attribute block `{...}` appears **after** the
bold caption text, which applies the styling to the `<p>` element (the
caption), not the `<table>` element.

### What Happened During HTML Generation

When Pandoc/knitr converted the markdown to HTML:

``` html
<!-- The table gets NO styling -->
<table>
  <thead>...</thead>
  <tbody>...</tbody>
</table>

<!-- Only the caption paragraph gets styled -->
<p style="width: 100%; margin-left: auto; margin-right: auto;">
  <strong>Table 3 | ...</strong> Divergence...
</p>
```

Result: The table retained default narrow width, while the caption had
100% width.

## Solution

Wrapped each table in a **Pandoc fenced div** (`::: ... :::`), which
properly contains the table and applies responsive styling.

### Fixed Markdown Structure

``` markdown
::: {.table-responsive style="width: 100%; overflow-x: auto;"}

|         |  q_0.01 |  q_0.05 |   q_0.1 |
|:--------|--------:|--------:|--------:|
| FOXJ2   | 0.04377 | 0.04633 | 0.04952 |
...

**Table 3 \| Pairwise Tsallis divergence estimates.** Divergence
(distance) between conditions...

:::
```

### Why This Works

1.  **Fenced div wrapper** (`:::`): Creates a
    `<div class="table-responsive">` container
2.  **width: 100%**: Table expands to full container width
3.  **overflow-x: auto**: Provides horizontal scrolling on narrow
    screens
4.  **Caption inside div**: Stays associated with the table visually

Generated HTML:

``` html
<div class="table-responsive" style="width: 100%; overflow-x: auto;">
  <table>
    <thead>...</thead>
    <tbody>...</tbody>
  </table>
  
  <p><strong>Table 3 | ...</strong> Divergence...</p>
</div>
```

## Files Modified

- `articles/TSENAT.md`
  - Table 3 (divergence estimates): lines ~847-870
  - Table 4 (effect sizes): lines ~886-910

## Changes Applied

| Table | Before | After |
|----|----|----|
| Both | Attributes on caption paragraph | Attributes on fenced div wrapper |
| Both | No scroll container | Responsive div with `overflow-x: auto` |

## Testing

After regenerating TSENAT.html with pkgdown:

``` bash
pkgdown::build_site()
```

✅ Tables now expand horizontally across full viewport width ✅
Mobile-responsive scrolling works on narrow screens ✅ Caption text
remains readable and properly associated

## Prevention

For future markdown tables in knitr documents:

**❌ Don’t do this:**

``` markdown
| col1 | col2 |
| ... | ... |

**Caption text** {style="width: 100%;"}
```

**✅ Do this:**

``` markdown
::: {style="width: 100%; overflow-x: auto;"}

| col1 | col2 |
| ... | ... |

**Caption text**

:::
```

## Related Issues

- Pandoc documentation: [Fenced
  Divs](https://pandoc.org/MANUAL.html#divs-and-spans)
- Bootstrap table CSS: Tables require parent container width to expand
- knitr/rmarkdown: CSS applied to wrong structural element leads to
  layout failures
