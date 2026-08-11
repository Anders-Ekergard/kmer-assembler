## name: test-every-function  
description: Use when writing or reviewing code with AI and you want it verified rather than just running. Makes the AI define success and rigorously test every function and its edge cases.


## Testing `inchworm_mvp.py`

### Functions to Test

The `inchworm_mvp.py` file contains the following functions that need testing:

- `most_common_kmers(kmers: list[str]) -> dict[str, int]`
- `kmers_in_contig(contig: str, k: int) -> set[str]`
- `expand_forward(kmers_count: dict[str, int], contig: str, k: int) -> str`
- `expand_backward(kmers_count: dict[str, int], contig: str, k: int) -> str`
- `expand(kmers_count: dict[str, int], k: int) -> list[str]`

### Test Cases

#### `most_common_kmers`

- **Success Criteria**: Returns a dictionary of k-mers and their counts, sorted by frequency.
- **Test Cases**:
  - **Ordinary Case**: Input: `['ATG', 'TGC', 'ATG', 'GCA']`, Expected: `{'ATG': 2, 'TGC': 1, 'GCA': 1}`
  - **Empty Input**: Input: `[]`, Expected: `{}`
  - **Single Input**: Input: `['ATG']`, Expected: `{'ATG': 1}`

#### `kmers_in_contig`

- **Success Criteria**: Returns all k-mers of length `k` in the contig.
- **Test Cases**:
  - **Ordinary Case**: Input: `contig = 'ATGCATGC'`, `k = 3`, Expected: `{'ATG', 'TGC', 'GCA', 'CAT', 'ATG'}`
  - **Edge Case**: Input: `contig = 'AT'`, `k = 3`, Expected: `set()`
  - **Single k-mer**: Input: `contig = 'ATG'`, `k = 3`, Expected: `{'ATG'}`

#### `expand_forward`

- **Success Criteria**: Expands the contig forward using the most common k-mers.
- **Test Cases**:
  - **Ordinary Case**: Input: `kmers_count = {'ATG': 2, 'TGC': 1, 'GCA': 1}`, `contig = 'ATG'`, `k = 3`, Expected: `'ATGC'` or similar
  - **No Expansion Possible**: Input: `kmers_count = {'ATG': 1}`, `contig = 'ATG'`, `k = 3`, Expected: `'ATG'`

#### `expand_backward`

- **Success Criteria**: Expands the contig backward using the most common k-mers.
- **Test Cases**:
  - **Ordinary Case**: Input: `kmers_count = {'ATG': 2, 'TGC': 1, 'GCA': 1}`, `contig = 'TGC'`, `k = 3`, Expected: `'ATGC'` or similar
  - **No Expansion Possible**: Input: `kmers_count = {'TGC': 1}`, `contig = 'TGC'`, `k = 3`, Expected: `'TGC'`

#### `expand`

- **Success Criteria**: Returns a list of contigs expanded in both directions.
- **Test Cases**:
  - **Ordinary Case**: Input: `kmers_count = {'ATG': 2, 'TGC': 1, 'GCA': 1}`, `k = 3`, Expected: `['ATGC']` or similar
  - **Empty Input**: Input: `kmers_count = {}`, `k = 3`, Expected: `[]`
