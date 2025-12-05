# rankdifferencetest 2025.12.4

- Major changes
  - `rdt()`
    - Complete rewrite with new features and breaking changes.
    - No longer depends on `coin::wilcoxsign_test()`, making it orders of magnitude faster.
  
- Additions
  - `srt()` for the Wilcoxon signed-rank test.
  - Pseudomedian and confidence interval of observed data.
    - `pmedian()`
    - `rdpmedian()`
  - Vector data interface
    - `rdt2()`
    - `srt2()`
    - `pmedian2()`
    - `rdpmedian2()`
  - Coerce objects to `data.frame` in typical `broom::tidy()` style.
    - `as.data.frame.rdt()`
    - `as.data.frame.srt()`
    - `as.data.frame.rdpmedian()`
    - `as.data.frame.pmedian()`
  - Coerce objects to tibble in typical `broom::tidy()` style.
    - `tidy.rdt()`
    - `tidy.srt()`
    - `tidy.rdpmedian()`
    - `tidy.pmedian()`
  - Brief overview of results in `htest` style.
    - `print.rdt()`
    - `print.srt()`
    - `print.pmedian()`
    - `print.rdpmedian()`
  
# rankdifferencetest 2025.4.21

- Fixed Rd `\link{}` targets missing package anchors.

# rankdifferencetest 2021-11-25

- First release of `rankdifferencetest`.
