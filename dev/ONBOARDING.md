# Claire's Onboarding Path

**Audience**: Claire — second-year PhD student, decent at R, new to package development

**Timeline**: 7 weeks (can be compressed or extended)

---

## Week 1-2: Understand the Data Structure

### Goals
- Understand what a `commaData` object is and what it contains
- Learn the accessor functions
- Get comfortable with the example dataset

### Tasks
1. Read `R/commaData_class.R` and `R/commaData_constructor.R` (don't worry about every detail)
2. In R console:
   ```r
   library(commaKit)
   data(comma_example_data)
   comma_example_data
   ```
3. Explore accessors:
   ```r
   methylation(comma_example_data)
   coverage(comma_example_data)
   modTypes(comma_example_data)
   modContexts(comma_example_data)
   siteInfo(comma_example_data)
   sampleInfo(comma_example_data)
   rowRanges(comma_example_data)  # genomic positions as GRanges
   ```
4. Read the "Getting Started" vignette:
   ```r
   vignette("getting-started", package = "commaKit")
   ```

### Questions to answer
- What's the difference between `mod_type` and `mod_context`?
- Why are genomic positions stored in `rowRanges()` instead of rowData?
- What does the `coverage()` matrix represent?

---

## Week 3-4: Understand the Analysis Pipeline

### Goals
- Walk through the core analysis functions
- Understand how methylation sites get annotated and tested

### Tasks
1. Run the full pipeline on example data:
   ```r
   # This is already in the vignette, but type it out yourself
   cdata <- comma_example_data

   # Annotate sites to genomic features
   cdata <- annotateSites(cdata, keep = "overlap")

   # Differential methylation testing
   cdata <- diffMethyl(cdata, formula = ~ condition, mod_type = "6mA")

   # Extract results
   res <- results(cdata)
   res_sig <- filterResults(cdata, padj = 0.05, delta_beta = 0.2)

   # Visualize
   plot_volcano(res_sig)
   plot_heatmap(res_sig, cdata, n_sites = 30L)
   ```

2. Try different methods:
   ```r
   cdata <- diffMethyl(cdata, formula = ~ condition, mod_type = "6mA", method = "methylkit")
   cdata <- diffMethyl(cdata, formula = ~ condition, mod_type = "6mA", method = "limma")
   ```

3. Compare results between methods

### Questions to answer
- What's the difference between the three `diffMethyl()` backends?
- Why might you choose one method over another?
- What does `annotateSites(keep = "overlap")` do vs. `keep = "proximity"`?

---

## Week 5-6: Package Development Basics

### Goals
- Learn roxygen2 documentation
- Learn testthat testing
- Understand the devtools workflow

### Tasks
1. Documentation:
   - Open `R/annotateSites.R`
   - Find the roxygen2 comments above the function
   - Understand what each tag does: `@param`, `@return`, `@export`, `@examples`
   - Run `devtools::document()` and see what changes in `man/`

2. Testing:
   - Open `tests/testthat/test-annotateSites.R`
   - Read through a few tests
   - Run `devtools::test()` and see the output
   - Try adding a simple test for a function you understand

3. Checks:
   - Run `devtools::check()` and read the output
   - Understand what it's checking

### Questions to answer
- Why do we write tests?
- What's the difference between `@export` and not exporting?
- What happens if you change a function but forget to run `document()`?

---

## Week 7+: Contributing

### Goals
- Make changes on a branch
- Write tests for new functionality
- Submit a PR

### Tasks
1. Create a feature branch:
   ```bash
   git checkout -b feature/my-first-change
   ```

2. Make a small change (e.g., add a new argument with a sensible default)

3. Add a test for it

4. Run `devtools::test()` and `devtools::check()`

5. Commit and push:
   ```bash
   git add -A
   git commit -m "Add new argument to function X"
   git push origin feature/my-first-change
   ```

6. Open a PR on GitHub

---

## Resources

- **commaBot**: Ask questions anytime — it's designed to help you learn
- **R Packages book**: <https://r-pkgs.org/> — excellent reference for package development
- **Bioconductor guidelines**: <https://contributions.bioconductor.org/>
- **testthat documentation**: <https://testthat.r-lib.org/>
- **roxygen2 documentation**: <https://roxygen2.r-lib.org/>

---

## Tips

- Don't rush. Understanding the structure is more important than speed.
- Ask commaBot questions — that's what it's for.
- If something is confusing, say so. The documentation probably needs improvement.
- Type code out yourself instead of copy-pasting. It helps memory.
