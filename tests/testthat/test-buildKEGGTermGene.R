## tests/testthat/test-buildKEGGTermGene.R
## Tests for buildKEGGTermGene() and the kegg_term2gene/kegg_term2name pathway
## in enrichMethylation().
##
## KEGGREST is NOT called for real -- all tests use local_mocked_bindings() to
## intercept KEGGREST::keggLink() and KEGGREST::keggList(), so no network
## access is required.

library(testthat)
library(comma)

# ── Shared mock data ──────────────────────────────────────────────────────────

# Mimics KEGGREST::keggLink("pathway", "eco") output
# names = gene IDs with organism prefix, values = pathway IDs with path: prefix
mock_links <- c(
    "eco:b0001" = "path:eco00010",
    "eco:b0002" = "path:eco00010",
    "eco:b0003" = "path:eco00020",
    "eco:b0004" = "path:eco00020",
    "eco:b0005" = "path:eco00030"
)

# Mimics KEGGREST::keggList("pathway", "eco") output
# names = pathway IDs with path: prefix, values = descriptions with organism suffix
mock_pathways <- c(
    "path:eco00010" = "Glycolysis / Gluconeogenesis - Escherichia coli K-12",
    "path:eco00020" = "Citrate cycle (TCA cycle) - Escherichia coli K-12",
    "path:eco00030" = "Pentose phosphate pathway - Escherichia coli K-12"
)

# ── buildKEGGTermGene() basic output ──────────────────────────────────────────

test_that("buildKEGGTermGene returns list with term2gene and term2name", {
    skip_if_not_installed("KEGGREST")
    local_mocked_bindings(
        keggLink = function(...) mock_links,
        keggList = function(...) mock_pathways,
        .package = "KEGGREST"
    )
    result <- buildKEGGTermGene("eco")
    expect_type(result, "list")
    expect_named(result, c("term2gene", "term2name"))
})

test_that("term2gene has columns 'term' and 'gene'", {
    skip_if_not_installed("KEGGREST")
    local_mocked_bindings(
        keggLink = function(...) mock_links,
        keggList = function(...) mock_pathways,
        .package = "KEGGREST"
    )
    result <- buildKEGGTermGene("eco")
    expect_s3_class(result$term2gene, "data.frame")
    expect_true(all(c("term", "gene") %in% colnames(result$term2gene)))
})

test_that("term2name has columns 'term' and 'name'", {
    skip_if_not_installed("KEGGREST")
    local_mocked_bindings(
        keggLink = function(...) mock_links,
        keggList = function(...) mock_pathways,
        .package = "KEGGREST"
    )
    result <- buildKEGGTermGene("eco")
    expect_s3_class(result$term2name, "data.frame")
    expect_true(all(c("term", "name") %in% colnames(result$term2name)))
})

test_that("strip_prefix = TRUE removes organism and path: prefixes", {
    skip_if_not_installed("KEGGREST")
    local_mocked_bindings(
        keggLink = function(...) mock_links,
        keggList = function(...) mock_pathways,
        .package = "KEGGREST"
    )
    result <- buildKEGGTermGene("eco", strip_prefix = TRUE)
    expect_false(any(grepl("^eco:", result$term2gene$gene)))
    expect_false(any(grepl("^path:", result$term2gene$term)))
    expect_false(any(grepl("^path:", result$term2name$term)))
})

test_that("strip_prefix = TRUE strips trailing organism qualifier from names", {
    skip_if_not_installed("KEGGREST")
    local_mocked_bindings(
        keggLink = function(...) mock_links,
        keggList = function(...) mock_pathways,
        .package = "KEGGREST"
    )
    result <- buildKEGGTermGene("eco", strip_prefix = TRUE)
    expect_false(any(grepl("Escherichia coli", result$term2name$name)))
    expect_true(any(grepl("Glycolysis", result$term2name$name)))
})

test_that("strip_prefix = FALSE retains all prefixes", {
    skip_if_not_installed("KEGGREST")
    local_mocked_bindings(
        keggLink = function(...) mock_links,
        keggList = function(...) mock_pathways,
        .package = "KEGGREST"
    )
    result <- buildKEGGTermGene("eco", strip_prefix = FALSE)
    expect_true(all(grepl("^eco:", result$term2gene$gene)))
    expect_true(all(grepl("^path:", result$term2gene$term)))
})

test_that("term2gene has correct number of rows", {
    skip_if_not_installed("KEGGREST")
    local_mocked_bindings(
        keggLink = function(...) mock_links,
        keggList = function(...) mock_pathways,
        .package = "KEGGREST"
    )
    result <- buildKEGGTermGene("eco")
    expect_equal(nrow(result$term2gene), length(mock_links))
})

test_that("term2name has correct number of rows", {
    skip_if_not_installed("KEGGREST")
    local_mocked_bindings(
        keggLink = function(...) mock_links,
        keggList = function(...) mock_pathways,
        .package = "KEGGREST"
    )
    result <- buildKEGGTermGene("eco")
    expect_equal(nrow(result$term2name), length(mock_pathways))
})

test_that("duplicate gene-pathway pairs are removed", {
    skip_if_not_installed("KEGGREST")
    dup_links <- c(mock_links, mock_links[1])  # add a duplicate
    local_mocked_bindings(
        keggLink = function(...) dup_links,
        keggList = function(...) mock_pathways,
        .package = "KEGGREST"
    )
    result <- buildKEGGTermGene("eco")
    expect_equal(nrow(result$term2gene), length(mock_links))
})

# ── Input validation ──────────────────────────────────────────────────────────

test_that("buildKEGGTermGene errors on non-character organism", {
    skip_if_not_installed("KEGGREST")
    expect_error(buildKEGGTermGene(123), "'organism' must be a non-empty character")
})

test_that("buildKEGGTermGene errors on empty-string organism", {
    skip_if_not_installed("KEGGREST")
    expect_error(buildKEGGTermGene(""), "'organism' must be a non-empty character")
})

test_that("buildKEGGTermGene errors on length-0 links", {
    skip_if_not_installed("KEGGREST")
    local_mocked_bindings(
        keggLink = function(...) character(0),
        keggList = function(...) mock_pathways,
        .package = "KEGGREST"
    )
    expect_error(buildKEGGTermGene("xxx"), "No KEGG pathway associations")
})

test_that("buildKEGGTermGene wraps keggLink errors with helpful message", {
    skip_if_not_installed("KEGGREST")
    local_mocked_bindings(
        keggLink = function(...) stop("HTTP 404"),
        keggList = function(...) mock_pathways,
        .package = "KEGGREST"
    )
    expect_error(buildKEGGTermGene("bad_org"), "KEGG API call failed")
})

test_that("buildKEGGTermGene warns when keggList fails and returns empty term2name", {
    skip_if_not_installed("KEGGREST")
    local_mocked_bindings(
        keggLink = function(...) mock_links,
        keggList = function(...) stop("server error"),
        .package = "KEGGREST"
    )
    result <- suppressWarnings(buildKEGGTermGene("eco"))
    expect_warning(buildKEGGTermGene("eco"), "Could not fetch KEGG pathway names")
    expect_equal(nrow(result$term2name), 0L)
})

# ── File caching ──────────────────────────────────────────────────────────────

test_that("buildKEGGTermGene saves and loads from RDS cache", {
    skip_if_not_installed("KEGGREST")
    tmp <- tempfile(fileext = ".rds")
    on.exit(unlink(tmp))

    local_mocked_bindings(
        keggLink = function(...) mock_links,
        keggList = function(...) mock_pathways,
        .package = "KEGGREST"
    )
    # First call: fetch and save
    r1 <- buildKEGGTermGene("eco", file = tmp)
    expect_true(file.exists(tmp))

    # Second call: load from cache (mock should NOT be called)
    called <- 0L
    local_mocked_bindings(
        keggLink = function(...) { called <<- called + 1L; mock_links },
        keggList = function(...) { called <<- called + 1L; mock_pathways },
        .package = "KEGGREST"
    )
    r2 <- buildKEGGTermGene("eco", file = tmp)
    expect_equal(called, 0L)
    expect_equal(r1$term2gene, r2$term2gene)
})

test_that("buildKEGGTermGene warns when cache file is older than 90 days", {
    skip_if_not_installed("KEGGREST")
    tmp <- tempfile(fileext = ".rds")
    on.exit(unlink(tmp))

    local_mocked_bindings(
        keggLink = function(...) mock_links,
        keggList = function(...) mock_pathways,
        .package = "KEGGREST"
    )
    buildKEGGTermGene("eco", file = tmp)

    # Backdate the file modification time by 100 days
    old_time <- Sys.time() - (100 * 24 * 3600)
    Sys.setFileTime(tmp, old_time)

    expect_warning(
        buildKEGGTermGene("eco", file = tmp),
        "days old"
    )
})

# ── Integration with enrichMethylation() ─────────────────────────────────────

test_that("enrichMethylation accepts kegg_term2gene and returns $kegg slot", {
    skip_if_not_installed("clusterProfiler")

    data(comma_example_data)
    dm  <- diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA")
    ann <- annotateSites(dm, annotation(comma_example_data), keep = "overlap")

    # Build a fake KEGG term2gene aligned with example data gene names
    fake_kegg_t2g <- data.frame(
        term = c("eco00010", "eco00010", "eco00020", "eco00020", "eco00030"),
        gene = c("geneA",    "geneB",    "geneC",    "geneD",    "geneE"),
        stringsAsFactors = FALSE
    )
    fake_kegg_t2n <- data.frame(
        term = c("eco00010", "eco00020", "eco00030"),
        name = c("Pathway Alpha", "Pathway Beta", "Pathway Gamma"),
        stringsAsFactors = FALSE
    )

    res <- enrichMethylation(
        ann,
        kegg_term2gene = fake_kegg_t2g,
        kegg_term2name = fake_kegg_t2n
    )

    expect_type(res, "list")
    expect_true("kegg" %in% names(res))
})

test_that("enrichMethylation kegg_term2gene result is in $kegg (not $go) slot", {
    skip_if_not_installed("clusterProfiler")

    data(comma_example_data)
    dm  <- diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA")
    ann <- annotateSites(dm, annotation(comma_example_data), keep = "overlap")

    fake_kegg_t2g <- data.frame(
        term = c("eco00010", "eco00010", "eco00020"),
        gene = c("geneA",    "geneB",    "geneC"),
        stringsAsFactors = FALSE
    )

    res <- enrichMethylation(ann, kegg_term2gene = fake_kegg_t2g, minGSSize = 1L)
    # $go should be NULL (no OrgDb/TERM2GENE), $kegg should be populated (even if NULL result from enricher)
    expect_null(res$go)
    # kegg slot exists in the result list (value may be NULL if no enrichment found)
    expect_true("kegg" %in% names(res))
})

test_that("enrichMethylation errors if kegg_term2gene missing required columns", {
    skip_if_not_installed("clusterProfiler")

    data(comma_example_data)
    dm  <- diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA")
    ann <- annotateSites(dm, annotation(comma_example_data), keep = "overlap")

    bad_t2g <- data.frame(pathway = "eco00010", gene_id = "geneA")

    expect_error(
        enrichMethylation(ann, kegg_term2gene = bad_t2g),
        "columns 'term' and 'gene'"
    )
})

test_that("enrichMethylation errors if kegg_term2name missing required columns", {
    skip_if_not_installed("clusterProfiler")

    data(comma_example_data)
    dm  <- diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA")
    ann <- annotateSites(dm, annotation(comma_example_data), keep = "overlap")

    fake_kegg_t2g <- data.frame(
        term = "eco00010", gene = "geneA",
        stringsAsFactors = FALSE
    )
    bad_t2n <- data.frame(pathway = "eco00010", description = "Glycolysis")

    expect_error(
        enrichMethylation(ann, kegg_term2gene = fake_kegg_t2g,
                          kegg_term2name = bad_t2n),
        "columns 'term' and 'name'"
    )
})

test_that("enrichMethylation error mentions kegg_term2gene in no-mapping message", {
    data(comma_example_data)
    dm  <- diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA")
    ann <- annotateSites(dm, annotation(comma_example_data), keep = "overlap")

    expect_error(
        enrichMethylation(ann),
        "kegg_term2gene"
    )
})

test_that("enrichMethylation kegg_term2gene works with method = 'gsea'", {
    skip_if_not_installed("clusterProfiler")

    data(comma_example_data)
    dm  <- diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA")
    ann <- annotateSites(dm, annotation(comma_example_data), keep = "overlap")

    fake_kegg_t2g <- data.frame(
        term = c("eco00010", "eco00010", "eco00020", "eco00020", "eco00030"),
        gene = c("geneA",    "geneB",    "geneC",    "geneD",    "geneE"),
        stringsAsFactors = FALSE
    )

    res <- enrichMethylation(
        ann,
        method         = "gsea",
        kegg_term2gene = fake_kegg_t2g,
        minGSSize      = 1L
    )
    expect_type(res, "list")
    expect_true("kegg" %in% names(res))
})

test_that("enrichMethylation kegg_term2gene works with method = c('ora','gsea')", {
    skip_if_not_installed("clusterProfiler")

    data(comma_example_data)
    dm  <- diffMethyl(comma_example_data, formula = ~ condition, mod_type = "6mA")
    ann <- annotateSites(dm, annotation(comma_example_data), keep = "overlap")

    fake_kegg_t2g <- data.frame(
        term = c("eco00010", "eco00010", "eco00020", "eco00020", "eco00030"),
        gene = c("geneA",    "geneB",    "geneC",    "geneD",    "geneE"),
        stringsAsFactors = FALSE
    )

    res <- enrichMethylation(
        ann,
        method         = c("ora", "gsea"),
        kegg_term2gene = fake_kegg_t2g,
        minGSSize      = 1L
    )
    expect_type(res, "list")
    expect_true("kegg" %in% names(res))
    expect_true(all(c("ora", "gsea") %in% names(res$kegg)))
})

# ── buildKEGGGeneIDMap() — shared mock data ───────────────────────────────────

# Mimics KEGGREST::keggConv("eco", "ncbi-geneid") output
# names = "ncbi-geneid:XXXXXX", values = "eco:bNNNN"
mock_conv <- c(
    "ncbi-geneid:945076" = "eco:b0344",  # lacZ
    "ncbi-geneid:945803" = "eco:b0343",  # lacY
    "ncbi-geneid:947498" = "eco:b0001",  # thrL
    "ncbi-geneid:945702" = "eco:b0002",  # thrA
    "ncbi-geneid:946367" = "eco:b0003"   # thrB
)

# Mimics a simple entrez2symbol data frame
mock_ent2sym <- data.frame(
    entrez_id = c("945076", "945803", "947498", "945702", "946367"),
    symbol    = c("lacZ",   "lacY",   "thrL",   "thrA",   "thrB"),
    stringsAsFactors = FALSE
)

# ── buildKEGGGeneIDMap() — entrez2symbol path ─────────────────────────────────

test_that("buildKEGGGeneIDMap entrez2symbol returns data.frame with symbol and kegg_id", {
    skip_if_not_installed("KEGGREST")
    local_mocked_bindings(
        keggConv = function(...) mock_conv,
        .package = "KEGGREST"
    )
    result <- buildKEGGGeneIDMap("eco", entrez2symbol = mock_ent2sym)
    expect_s3_class(result, "data.frame")
    expect_true(all(c("symbol", "kegg_id") %in% colnames(result)))
})

test_that("buildKEGGGeneIDMap entrez2symbol kegg_id has no organism prefix", {
    skip_if_not_installed("KEGGREST")
    local_mocked_bindings(
        keggConv = function(...) mock_conv,
        .package = "KEGGREST"
    )
    result <- buildKEGGGeneIDMap("eco", entrez2symbol = mock_ent2sym)
    expect_false(any(grepl("^eco:", result$kegg_id)))
    expect_false(any(grepl("^ncbi-geneid:", result$kegg_id)))
})

test_that("buildKEGGGeneIDMap entrez2symbol join is correct", {
    skip_if_not_installed("KEGGREST")
    local_mocked_bindings(
        keggConv = function(...) mock_conv,
        .package = "KEGGREST"
    )
    result <- buildKEGGGeneIDMap("eco", entrez2symbol = mock_ent2sym)
    lacz_row <- result[result$symbol == "lacZ", ]
    expect_equal(nrow(lacz_row), 1L)
    expect_equal(lacz_row$kegg_id, "b0344")
})

test_that("buildKEGGGeneIDMap entrez2symbol excludes genes absent from keggConv", {
    skip_if_not_installed("KEGGREST")
    extra_sym <- rbind(
        mock_ent2sym,
        data.frame(entrez_id = "999999", symbol = "fakeGene",
                   stringsAsFactors = FALSE)
    )
    local_mocked_bindings(
        keggConv = function(...) mock_conv,
        .package = "KEGGREST"
    )
    result <- buildKEGGGeneIDMap("eco", entrez2symbol = extra_sym)
    expect_false("fakeGene" %in% result$symbol)
    expect_equal(nrow(result), nrow(mock_ent2sym))
})

test_that("buildKEGGGeneIDMap errors when entrez2symbol missing required columns", {
    skip_if_not_installed("KEGGREST")
    bad <- data.frame(gene_id = "945076", name = "lacZ")
    expect_error(
        buildKEGGGeneIDMap("eco", entrez2symbol = bad),
        "'entrez_id' and 'symbol'"
    )
})

test_that("buildKEGGGeneIDMap errors when neither OrgDb nor entrez2symbol provided", {
    skip_if_not_installed("KEGGREST")
    expect_error(
        buildKEGGGeneIDMap("eco"),
        "Provide at least one of"
    )
})

test_that("buildKEGGGeneIDMap errors on bad organism", {
    skip_if_not_installed("KEGGREST")
    expect_error(buildKEGGGeneIDMap(123, entrez2symbol = mock_ent2sym),
                 "'organism' must be a non-empty character")
    expect_error(buildKEGGGeneIDMap("",  entrez2symbol = mock_ent2sym),
                 "'organism' must be a non-empty character")
})

test_that("buildKEGGGeneIDMap wraps keggConv errors with helpful message", {
    skip_if_not_installed("KEGGREST")
    local_mocked_bindings(
        keggConv = function(...) stop("HTTP 404"),
        .package = "KEGGREST"
    )
    expect_error(
        buildKEGGGeneIDMap("bad_org", entrez2symbol = mock_ent2sym),
        "KEGG API call failed"
    )
})

test_that("buildKEGGGeneIDMap warns when no symbol matches found", {
    skip_if_not_installed("KEGGREST")
    no_match <- data.frame(entrez_id = "000000", symbol = "X",
                           stringsAsFactors = FALSE)
    local_mocked_bindings(
        keggConv = function(...) mock_conv,
        .package = "KEGGREST"
    )
    expect_warning(
        buildKEGGGeneIDMap("eco", entrez2symbol = no_match),
        "No symbol matches found"
    )
})

# ── buildKEGGGeneIDMap() — OrgDb path ─────────────────────────────────────────

# Build a minimal mock OrgDb-like environment
mock_orgdb <- new.env()

test_that("buildKEGGGeneIDMap OrgDb path returns correct symbol-kegg_id pairs", {
    skip_if_not_installed("KEGGREST")
    skip_if_not_installed("AnnotationDbi")
    local_mocked_bindings(
        keggConv = function(...) mock_conv,
        .package = "KEGGREST"
    )
    local_mocked_bindings(
        columns = function(x, ...) c("ENTREZID", "SYMBOL"),
        keys    = function(x, keytype, ...) mock_ent2sym$entrez_id,
        select  = function(x, keys, columns, keytype, ...) {
            data.frame(
                ENTREZID = mock_ent2sym$entrez_id,
                SYMBOL   = mock_ent2sym$symbol,
                stringsAsFactors = FALSE
            )
        },
        .package = "AnnotationDbi"
    )
    result <- buildKEGGGeneIDMap("eco", OrgDb = mock_orgdb)
    expect_s3_class(result, "data.frame")
    expect_true(all(c("symbol", "kegg_id") %in% colnames(result)))
    lacz_row <- result[result$symbol == "lacZ", ]
    expect_equal(lacz_row$kegg_id, "b0344")
})

test_that("buildKEGGGeneIDMap OrgDb path errors when id_col not in OrgDb", {
    skip_if_not_installed("KEGGREST")
    skip_if_not_installed("AnnotationDbi")
    local_mocked_bindings(
        keggConv = function(...) mock_conv,
        .package = "KEGGREST"
    )
    local_mocked_bindings(
        columns = function(x, ...) c("SYMBOL"),  # ENTREZID missing
        .package = "AnnotationDbi"
    )
    expect_error(
        buildKEGGGeneIDMap("eco", OrgDb = mock_orgdb, id_col = "ENTREZID"),
        "not found in OrgDb"
    )
})

# ── buildKEGGGeneIDMap() — caching ────────────────────────────────────────────

test_that("buildKEGGGeneIDMap saves and loads from RDS cache", {
    skip_if_not_installed("KEGGREST")
    tmp <- tempfile(fileext = ".rds")
    on.exit(unlink(tmp))
    local_mocked_bindings(
        keggConv = function(...) mock_conv,
        .package = "KEGGREST"
    )
    # First call: fetch and save
    r1 <- buildKEGGGeneIDMap("eco", entrez2symbol = mock_ent2sym, file = tmp)
    expect_true(file.exists(tmp))

    # Second call: load from cache — keggConv must NOT be called
    called <- 0L
    local_mocked_bindings(
        keggConv = function(...) { called <<- called + 1L; mock_conv },
        .package = "KEGGREST"
    )
    r2 <- buildKEGGGeneIDMap("eco", entrez2symbol = mock_ent2sym, file = tmp)
    expect_equal(called, 0L)
    expect_equal(r1, r2)
})

test_that("buildKEGGGeneIDMap warns when cache older than 90 days", {
    skip_if_not_installed("KEGGREST")
    tmp <- tempfile(fileext = ".rds")
    on.exit(unlink(tmp))
    local_mocked_bindings(
        keggConv = function(...) mock_conv,
        .package = "KEGGREST"
    )
    buildKEGGGeneIDMap("eco", entrez2symbol = mock_ent2sym, file = tmp)
    Sys.setFileTime(tmp, Sys.time() - (100 * 24 * 3600))
    expect_warning(
        buildKEGGGeneIDMap("eco", entrez2symbol = mock_ent2sym, file = tmp),
        "days old"
    )
})

test_that("buildKEGGGeneIDMap with file=NULL does not write a file", {
    skip_if_not_installed("KEGGREST")
    local_mocked_bindings(
        keggConv = function(...) mock_conv,
        .package = "KEGGREST"
    )
    old_files <- list.files(tempdir(), pattern = "\\.rds$")
    buildKEGGGeneIDMap("eco", entrez2symbol = mock_ent2sym, file = NULL)
    new_files <- list.files(tempdir(), pattern = "\\.rds$")
    expect_equal(length(new_files), length(old_files))
})

# ── buildKEGGTermGene(id_map = ...) integration ───────────────────────────────

# id_map covering b0001-b0005 used by mock_links
mock_id_map <- data.frame(
    symbol  = c("thrL", "thrA", "aceE", "aceF", "lpdA"),
    kegg_id = c("b0001", "b0002", "b0003", "b0004", "b0005"),
    stringsAsFactors = FALSE
)

test_that("buildKEGGTermGene with id_map translates gene column to symbols", {
    skip_if_not_installed("KEGGREST")
    local_mocked_bindings(
        keggLink = function(...) mock_links,
        keggList = function(...) mock_pathways,
        .package = "KEGGREST"
    )
    result <- buildKEGGTermGene("eco", id_map = mock_id_map)
    expect_true(all(result$term2gene$gene %in% mock_id_map$symbol))
    expect_false(any(grepl("^b[0-9]", result$term2gene$gene)))
})

test_that("buildKEGGTermGene with id_map preserves unmatched KEGG IDs", {
    skip_if_not_installed("KEGGREST")
    partial_map <- mock_id_map[1:3, ]  # only b0001-b0003
    local_mocked_bindings(
        keggLink = function(...) mock_links,
        keggList = function(...) mock_pathways,
        .package = "KEGGREST"
    )
    result <- buildKEGGTermGene("eco", id_map = partial_map)
    genes <- result$term2gene$gene
    # b0001-b0003 should be translated; b0004 and b0005 preserved as-is
    expect_true(all(c("thrL", "thrA", "aceE") %in% genes))
    expect_true(all(c("b0004", "b0005") %in% genes))
    # No NAs
    expect_false(any(is.na(genes)))
})

test_that("buildKEGGTermGene without id_map is unchanged from original behavior", {
    skip_if_not_installed("KEGGREST")
    local_mocked_bindings(
        keggLink = function(...) mock_links,
        keggList = function(...) mock_pathways,
        .package = "KEGGREST"
    )
    result <- buildKEGGTermGene("eco")
    # Default behavior: b-numbers in gene column
    expect_true(all(grepl("^b[0-9]", result$term2gene$gene)))
})

test_that("buildKEGGTermGene with id_map = NULL is identical to default", {
    skip_if_not_installed("KEGGREST")
    local_mocked_bindings(
        keggLink = function(...) mock_links,
        keggList = function(...) mock_pathways,
        .package = "KEGGREST"
    )
    r_default <- buildKEGGTermGene("eco")
    r_null    <- buildKEGGTermGene("eco", id_map = NULL)
    expect_equal(r_default, r_null)
})

test_that("buildKEGGTermGene errors on malformed id_map", {
    skip_if_not_installed("KEGGREST")
    bad_map <- data.frame(gene_name = "thrL", identifier = "b0001")
    local_mocked_bindings(
        keggLink = function(...) mock_links,
        keggList = function(...) mock_pathways,
        .package = "KEGGREST"
    )
    expect_error(
        buildKEGGTermGene("eco", id_map = bad_map),
        "'symbol' and 'kegg_id'"
    )
})
