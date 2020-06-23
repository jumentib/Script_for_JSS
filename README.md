# Sparse Latent Factor Regression Models for Association Studies in R

This project contains all the R scripts used to conduct simulation studies and real data studies presented in the article: Sparse Latent Factor Regression Models for Association Studies in R.

## Empirical simulation experiments

This section reproduces the experiments described in the section "Empirical simulation experiments" of the manuscript.

Empirical simulations are based on empirical genotype data for the model plant *Arabidopsis thaliana* from (Atwell et al. 2010). An artificial dependent variable, or phenotype, is simulated from five causal quantitative trait loci (causal loci).  

## Generative model simulations

This section reproduces the experiments described in section "Generative simulation experiments" of the manuscript. Generative simulations are based on the generative model of LFMM, explained in the main text equation (1).

## Genome-wide association study

We carried out genome-wide association studies (GWAS) between *Arabidopsis thaliana* genotypes and flowering time phenotypes (FT16). GWAS were performed with 6 different methods: sparse lfmm, bslmm, lasso, ridge lfmm, cate and sva.

## Data

This directory contains the genotype data of *Arabidopsis thaliana*.
Which contains n = 162 European accessions and p = 53,859 SNPs from the fifth chromosome of the Arabidopsis thaliana genome to investigate associations with a flowering time phenotype (FT16-TO: 0000344, (Atwell et al. 2010)). FT16 corresponds to the number of days required for an individual plant to reach the flowering stage.

## Gemma

This file contains an executable from the Gemma program which includes the Bayesian Sparse Linear Mixed Model (BSLMM) from Zhou et al (2013).

## Function

This script contains a set of functions for carrying out simulation studies. Functions for wrapping the Gemma program. Functions to calculate the power of the different statistical methods tested during the simulation studies.

