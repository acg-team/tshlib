# Tree-Search-Library Readme

The Tree-Search-Library (TSH-LIB) is a new implementation of binary tree topological searching engines

## Compiling

## Installing

## Integrating into an existing project

## Workflow

Given an utree as input, the algorithm performs the following steps:

1. Per each node in the tree (either internal or terminal), find 

## Algorithms

The library provides an automatic heuristic to sample the tree-space given a scoring function, an initial tree topology and an alignment.

 1. automatic tree scaling strategy
 2. swarm particle-like strategy

## Methods

## For developers

### Naming conventions

Class methods are named according to the following convention:
 1. public methods: verbObject (i.e. printNodes, listMoves)
 2. private methods: _verbObject (i.e. _testTopology)

Variables are named with a self explaining indication of the variable contents.

### Data structures

### Extending the library

### Benchmark

The library can be benchmarked using the utility tshlib-benchmark included in the directory ./benchmarks.
