# BUG
Bioinformatic utilities in Go

## Allelespecificexp

An oversimplified method to find allele specific expression in a normalized expression table and pedigree information

Usage: allelespecificexp -exp table.tsv -group "M M M F F F 1 1 1 2 2 2 3 3 3"

M marks the column of the mother RNA expression.

F marks the column of the father RNA expression.

Numbers marks the offsprings.

Replicates should have similar group ids.

## Fastagrep

Using regular expression to find fasta records

## Fastasubsample

Create a smaller subset of a fasta file

## hmm

Simple HMM framework. Right now it is just a test framework
