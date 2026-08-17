# PhyloDist

**PhyloDist** is a Python library for comparing phylogenetic trees using classical and information-theoretic tree distance and similarity measures.

The library is designed to:

- Reproduce results from **TreeDist (R)**
- Support both **rooted** and **unrooted** tree comparisons
- Compare rooted trees as unrooted trees when requested, by ignoring the root direction
- Provide both **distance** and **similarity** measures
- Support **normalised scores in [0,1]**
- Support optional taxon-name mapping before comparison
- Provide a clean and extensible Python interface
- Enable validation against published phylogenetic metrics

---

## Features

Supported comparison methods:

| Method | Description | Return Type | TreeDist (R) Equivalent |
|--------|-------------|-------------|--------------------------|
| `rf` | Robinson–Foulds distance | Distance | `RobinsonFoulds(tree1, tree2, normalise = FALSE)` |
| `jrf` | Jaccard Robinson–Foulds distance | Distance | `JaccardRobinsonFoulds(tree1, tree2, normalise = FALSE)` |
| `icrf` | Information-corrected RF | Similarity or TreeDist-equivalent distance | `InfoRobinsonFoulds(tree1, tree2, normalise = FALSE)` |
| `nye` | Nye similarity | Similarity | `NyeSimilarity(tree1, tree2, normalise = FALSE)` |
| `msi` | Matching Split Information | Similarity | `MatchingSplitInfo(tree1, tree2, normalise = FALSE)` |
| `msid` | Matching Split Information Distance | Distance | `MatchingSplitInfoDistance(tree1, tree2, normalise = FALSE)` |
| `spi` | Shared Phylogenetic Information | Similarity | `SharedPhylogeneticInfo(tree1, tree2, normalise = FALSE)` |
| `pid` | Phylogenetic Information Distance | Distance | `PhylogeneticInfoDistance(tree1, tree2, normalise = FALSE)` |
| `mci` | Mutual Clustering Information | Similarity | `MutualClusteringInfo(tree1, tree2, normalise = FALSE)` |
| `cid` | Clustering Information Distance | Distance | `ClusteringInfoDist(tree1, tree2, normalise = FALSE)` |
| `frobenius` | Frobenius distance between leaf path-length matrices | Distance | — |

Additional support:

- Plain Newick strings or tree files
- NEXUS-style tree blocks
- Files containing multiple Newick trees; one complete tree is selected before parsing
- ETE4 parsing with fallback topology-only parsing for difficult Newick annotations
- Optional taxon mapping, useful when gene-tree leaves and species-tree leaves use different naming conventions

### Branch-length distance

Frobenius distance compares the patristic (leaf-to-leaf path-length) matrices,
so it captures both topology and branch-length differences. Each unordered
leaf pair is counted once. The trees must have exactly the same uniquely named
leaves and must retain branch lengths through ETE4 parsing.

```python
comparison = PhyloDist.from_newick(tree1, tree2, parser=1)

raw = comparison.compare(method="frobenius").score
relative = comparison.compare(
    method="frobenius",
    frobenius_type="relative",
).score
```

The relative value is `raw / ||D2||`, so tree 2 is the reference and reversing
the trees can change the result. This option is separate from `normalise`
because a relative Frobenius distance is not restricted to `[0, 1]`.

---

## Installation

### Method 1 - Install directly from GitHub

```bash
pip install git+https://github.com/OrthoFinder/PhyloDist.git
```

### Method 2 - Clone the repository

```bash
git clone https://github.com/OrthoFinder/PhyloDist.git
cd PhyloDist
pip install .
```

### Dependencies

- Python ≥ 3.11
- numpy
- scipy
- ete4

---

## Basic Python Usage

```python
from phylodist import PhyloDist, Tree

# Trees may be Newick strings or file paths.
tree1 = "tree1.txt"
tree2 = "tree2.txt"

# Or:
# tree1 = "((A,B),(C,D));"
# tree2 = "((A,C),(B,D));"

t1 = Tree(
    tree1,
    is_rooted=False,
    parser=1,
)

t2 = Tree(
    tree2,
    is_rooted=False,
    parser=1,
)

grf = PhyloDist(t1, t2)

result = grf.compare(
    method="mci",
    rooted=False,      # False = unrooted split comparison
    normalise=True,    # Scale result to [0,1]
    return_type="similarity",
)

print(result.score)
```

For a distance-style example, use the distance alias `cid`:

```python
result = grf.compare(
    method="cid",
    rooted=False,
    normalise=True,
    return_type="distance",
)
```

---

## Rooted and Unrooted Comparisons

PhyloDist separates two concepts:

1. **Tree rootedness metadata**: whether each input tree is marked as rooted or unrooted.
2. **Comparison mode**: whether PhyloDist compares rooted clades or unrooted splits.

This matters because a rooted tree can be compared as an unrooted tree by ignoring the root direction, but an unrooted tree cannot safely be compared as rooted without a biologically meaningful root.

### Comparison modes

| Mode | Meaning | Recommended use |
|------|---------|-----------------|
| `rooted` | Compare rooted clades | Only when both trees are rooted using comparable biological rooting |
| `unrooted` | Compare unrooted splits | Recommended for most gene-tree comparisons |
| `auto` | Use rooted only if both trees are explicitly rooted; otherwise use unrooted | Useful as a general CLI default |

For gene/species tree benchmarking, such as rooted species tree versus unrooted inferred gene trees, use unrooted comparison.

---

## Taxon Mapping

PhyloDist can optionally map leaf names before comparison. This is useful when two trees use different naming conventions for the same taxa.

For example, a species tree may use full species names:

```text
Saccharomyces_cerevisiae
```

while a gene tree may use sequence names or species-code prefixes:

```text
Scer_YAL001C
```

With a taxon mapper, `Scer_YAL001C` can be mapped to:

```text
Saccharomyces_cerevisiae
```

before tree comparison.

### Supported mapper formats

PhyloDist accepts either whitespace/TSV format:

```text
pattern replacement
Scer    Saccharomyces_cerevisiae
Smik    Saccharomyces_mikatae
Skud    Saccharomyces_kudriavzevii
```

or colon format:

```text
Scer: Saccharomyces_cerevisiae
Smik: Saccharomyces_mikatae
Skud: Saccharomyces_kudriavzevii
```

The header is optional. These are also valid:

```text
Scer    Saccharomyces_cerevisiae
Smik    Saccharomyces_mikatae
Skud    Saccharomyces_kudriavzevii
```

```text
Scer: Saccharomyces_cerevisiae
Smik: Saccharomyces_mikatae
Skud: Saccharomyces_kudriavzevii
```

### Comments

Empty lines are ignored. Lines starting with `#` or `//` are ignored. Inline comments after `#` or `//` are also ignored:

```text
# Yeast taxon mapper
// Short code to full species name

Scer: Saccharomyces_cerevisiae  # budding yeast
Smik: Saccharomyces_mikatae     // mapped from gene-tree prefix
```

Optional headers such as the following are ignored:

```text
pattern replacement
from to
source target
old new
input output
```

### Matching behaviour

For each leaf name, PhyloDist tries mapping rules in this order:

1. exact full-name match;
2. prefix before the first underscore;
3. regular-expression full match;
4. unchanged name if no rule matches.

For example, with this mapper:

```text
Scer: Saccharomyces_cerevisiae
```

both of these leaves map to `Saccharomyces_cerevisiae`:

```text
Scer
Scer_YAL001C
```

### Python usage with taxon mapping

The `taxon_mapper` argument can be:

- `None`
- a mapper file path
- a callable function
- a list or tuple of `(pattern, replacement)` rules

PhyloDist resolves mapper files internally, so downstream tools can pass the mapper file path directly:

```python
from phylodist import PhyloDist, Tree

t1 = Tree(
    "species_tree.tree",
    is_rooted=True,
    parser=None,
    taxon_mapper="yeast_taxon_mapper.tsv",
)

t2 = Tree(
    "gene_tree.tree",
    is_rooted=False,
    parser=None,
    taxon_mapper="yeast_taxon_mapper.tsv",
)

grf = PhyloDist(t1, t2)

result = grf.compare(
    method="cid",
    rooted=False,
    normalise=True,
    return_type="distance",
)

print(result.score)
```

You can also pass mapping rules directly:

```python
rules = [
    ("Scer", "Saccharomyces_cerevisiae"),
    ("Smik", "Saccharomyces_mikatae"),
]

t = Tree(
    "gene_tree.tree",
    taxon_mapper=rules,
)
```

or use a custom callable:

```python
def mapper(name):
    return name.split("_", 1)[0]

t = Tree(
    "gene_tree.tree",
    taxon_mapper=mapper,
)
```

Taxon mapping is applied while trees are loaded, before rooted clades or unrooted splits are extracted.

---

## Command Line Usage

PhyloDist can be run in several ways depending on how it is installed.

### Method 1 - Installed command line tool

```bash
phylodist tree1.txt tree2.txt \
    --method cid \
    --tree1-unrooted \
    --tree2-unrooted \
    --compare-mode unrooted \
    --normalise
```

### Method 2 - Run using the standalone script

```bash
python phylodist.py tree1.txt tree2.txt \
    --method cid \
    --tree1-unrooted \
    --tree2-unrooted \
    --compare-mode unrooted \
    --normalise
```

### Example with taxon mapping

```bash
phylodist species_tree.tree gene_tree.tree \
    --method cid \
    --tree1-rooted \
    --tree2-unrooted \
    --compare-mode unrooted \
    --taxon-mapper yeast_taxon_mapper.tsv \
    --normalise
```

### Example using auto mode

```bash
phylodist tree1.txt tree2.txt \
    --method rf \
    --compare-mode auto
```

In `auto` mode, PhyloDist uses rooted comparison only when both trees are explicitly marked rooted. Otherwise, it uses unrooted comparison.

---

## Command Line Options

### Required Arguments

| Argument | Description |
|----------|-------------|
| `tree1` | First tree file or Newick string |
| `tree2` | Second tree file or Newick string |

### General Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--method` | string | `rf` | Comparison method |
| `--normalise` | flag | `False` | Return normalised score in `[0,1]` |
| `--return-type` | string | `distance` | Return `distance` or `similarity` |
| `--k` | float | `1.0` | Exponent parameter for JRF |
| `--frobenius-type` | string | `raw` | Return `raw` or tree-2-relative Frobenius distance |
| `--parser` | int | `1` | ETE parser format |
| `--taxon-mapper` | path | `None` | Optional taxon mapper file |

### Rootedness and Comparison Mode Options

| Option | Description |
|--------|-------------|
| `--tree1-rooted` | Mark tree1 as rooted |
| `--tree1-unrooted` | Mark tree1 as unrooted |
| `--tree2-rooted` | Mark tree2 as rooted |
| `--tree2-unrooted` | Mark tree2 as unrooted |
| `--compare-mode auto` | Use rooted comparison only if both trees are explicitly rooted; otherwise unrooted |
| `--compare-mode rooted` | Force rooted clade comparison |
| `--compare-mode unrooted` | Force unrooted split comparison |

### Output Options

| Option | Description |
|--------|-------------|
| `--verbose` | Print diagnostic information |
| `--print-matching` | Display matched split/clade pairs |
| `--json` | Output full result as JSON |

### Available Methods

| Method | Description | Output |
|--------|-------------|--------|
| `rf` | Robinson–Foulds | Distance |
| `jrf` | Jaccard Robinson–Foulds | Distance |
| `icrf` | Information-corrected RF | Similarity or Distance |
| `nye` | Nye similarity | Similarity |
| `msi` | Matching Split Information | Similarity |
| `msid` | Matching Split Information Distance | Distance |
| `spi` | Shared Phylogenetic Information | Similarity |
| `pid` | Phylogenetic Information Distance | Distance |
| `mci` | Mutual Clustering Information | Similarity |
| `cid` | Clustering Information Distance | Distance |

### Get Help

```bash
phylodist --help
```

---

## Scalability Notes

Unlike **TreeDist (R)**, which is often used on trees of approximately **≤1000 leaves** in practice, **PhyloDist** does not impose a fixed upper limit on tree size.

However, many tree comparison methods implemented here rely on bipartite matching algorithms such as the Hungarian algorithm, whose computational complexity grows rapidly with tree size.

Typical complexity:

| Component | Complexity |
|----------|-------------|
| Split comparison | O(n²) |
| Matching, Hungarian algorithm | O(n³) |
| Memory usage | O(n²) |

Where:

- `n` = number of splits, roughly proportional to number of leaves

### Practical Recommendations

For best performance:

- ≤ **500 leaves** — fast
- **500–2000 leaves** — moderate runtime
- ≥ **2000 leaves** — potentially slow
- ≥ **5000 leaves** — large memory usage expected

Performance depends strongly on:

- number of taxa
- tree balance
- chosen metric
- available system memory

Large trees are supported, but runtime may increase substantially. Users working with large phylogenies are encouraged to test on smaller subsets first, monitor memory usage, and consider parallel workflows for batch comparisons.

---

## Important Interpretation Notes

Many tree comparison metrics implemented here are information-based. For these metrics:

- Similarity values increase as trees become more similar.
- Distance values increase as trees become more different.
- Raw information-based scores depend on tree size.
- These values should not be compared across datasets without normalization.
- Use `normalise=True` when comparing trees of different sizes.

Metrics especially affected include:

- `icrf`
- `msi` / `msid`
- `spi` / `pid`
- `mci` / `cid`

Rooted and unrooted comparisons are not interchangeable. Rooted mode compares rooted clades, while unrooted mode compares unrooted splits. Use rooted mode only when the input trees have comparable biological rooting.

---

## Reference

- Martin R Smith, Information theoretic generalized Robinson–Foulds metrics for comparing phylogenetic trees, Bioinformatics, Volume 36, Issue 20, October 2020, Pages 5007–5013, [![DOI](https://img.shields.io/badge/DOI-10.1093%2Fbioinformatics%2Fbtaa614-blue)](https://doi.org/10.1093/bioinformatics/btaa614)
