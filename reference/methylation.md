# Accessor for the methylation (beta value) matrix

Retrieves the sites × samples matrix of methylation beta values from a
[`commaData`](https://carl-stone.github.io/commaKit/reference/commaData.md)
object. Values are in the range 0–1 (proportion of reads called
methylated). Sites below the `min_coverage` threshold set at object
creation are `NA`.

## Usage

``` r
methylation(object)

# S4 method for class 'commaData'
methylation(object)

# S4 method for class 'ANY'
methylation(object)
```

## Arguments

- object:

  A `commaData` object.

## Value

A numeric matrix with rows corresponding to methylation sites and
columns corresponding to samples. Rownames are site keys
(`"chrom:position:strand:mod_type"`); column names are sample names.

## See also

[`siteCoverage`](https://carl-stone.github.io/commaKit/reference/siteCoverage.md),
[`siteInfo`](https://carl-stone.github.io/commaKit/reference/siteInfo.md),
[`sampleInfo`](https://carl-stone.github.io/commaKit/reference/sampleInfo.md)

## Examples

``` r
data(comma_example_data)
m <- methylation(comma_example_data)
dim(m)
#> [1] 588   6
head(m)
#>         ctrl_1    ctrl_2    ctrl_3   treat_1   treat_2   treat_3
#> [1,] 0.7575443 0.8107216 0.8870622 0.9077817 0.8509780 0.9272379
#> [2,] 0.9738114 0.9121264 0.8963508 0.9592666 0.8734184 0.7990160
#> [3,] 0.8678683 0.8581721 0.9193607 0.9801922 0.9301237 0.9180024
#> [4,] 0.9569443 0.9232934 0.8870824 0.8679993 0.9137143 0.9810536
#> [5,] 0.9165347 0.8104115 0.9441796 0.9399405 0.9162951 0.6952466
#> [6,] 0.8337892 0.9139209 0.9433962 0.4032709 0.3884889 0.3921812
```
