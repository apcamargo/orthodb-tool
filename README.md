# `orthodb-tool`

`orthodb-tool` is a command-line tool for querying [OrthoDB](https://www.orthodb.org/). It can export orthologous group records as TSV, summarize OG counts across descendant taxa, and download FASTA sequences for a single OG.

## Installation

You can execute `orthodb-tool` using [Pixi](https://pixi.sh/):

```bash
git clone https://github.com/apcamargo/orthodb-tool.git
cd orthodb-tool
pixi install
```

## Usage

```bash
pixi run orthodb-tool --help
pixi run orthodb-tool <subcommand> --help
```

`orthodb-tool` writes the output to `stdout` by default. Use `--output` to write to a file.

## Subcommands

### `counts`

Count OGs for descendant taxa under a parent taxid at a target rank. The TSV output includes both `n_genomes` and `n_og`.

```bash
pixi run orthodb-tool counts --parent-taxid 1236 --target-rank family --output counts.tsv
```

```
 taxid   taxon_name           n_genomes   n_og
────────────────────────────────────────────────
 641     Vibrionaceae         204         21702
 32033   Xanthomonadaceae     200         18131
 543     Enterobacteriaceae   175         17539
 72275   Alteromonadaceae     78          14949
 468     Moraxellaceae        175         14690
```

### `records`

Export orthologous group records as TSV. Optional flags can add descendant taxids, taxonomy rank, annotation ID columns, and parent OG IDs.

```bash
pixi run orthodb-tool records --taxids 32033 --min-genes 10 --universal 0.8 --output records.tsv
```

```
 og_id          taxid   taxon_name         protein_count   description
───────────────────────────────────────────────────────────────────────────────────────────────────
 10023at32033   32033   Xanthomonadaceae   169             AAA domain
 1004at32033    32033   Xanthomonadaceae   196             polyketide cyclase
 1009at32033    32033   Xanthomonadaceae   203             tetratricopeptide repeat protein
 1011at32033    32033   Xanthomonadaceae   199             thymidylate synthase
 1012at32033    32033   Xanthomonadaceae   201             ATP synthase, F1 complex, gamma subunit
```

### `download-fasta`

Download the sequences of proteins in an OG in FASTA format. By default, the tool downloads protein sequences, but you can use `--cds` to download CDS sequences instead.

```bash
pixi run orthodb-tool download-fasta 1037at32033 --output sequences.fasta
```
