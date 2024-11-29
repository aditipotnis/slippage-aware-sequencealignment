import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency

# Load the phenotypes
phenotypes = pd.read_csv('/Users/aditipotnis/Downloads/gwas_phenotypes.txt', sep='\t', header=None, names=['individual', 'phenotype'])

# Load the VCF file
with open('/Users/aditipotnis/Downloads/gwas_population.vcf', 'r') as f:
    lines = [line.strip() for line in f if not line.startswith('#')]

# Process VCF data (assuming each SNP is encoded in "0|0" format)
genotypes = []
snp_ids = []
for line in lines:
    parts = line.split('\t')
    snp_ids.append(parts[2])  # Assuming the 3rd column contains SNP IDs
    # Split each genotype like "0|0", "0|1", or "1|1"
    genotypes.append([sum(map(int, allele.split('|'))) if allele != '.' else np.nan for allele in parts[6:]])

genotypes = np.array(genotypes).T  # Transpose to have individuals as rows and SNPs as columns

# Chi-squared test for each SNP
p_values = []
odds_ratios_het = []
odds_ratios_homo_alt = []
num_snps = genotypes.shape[1]

for snp in range(num_snps):
    contingency_table = pd.crosstab(genotypes[:, snp], phenotypes['phenotype'])
    chi2, p, dof, expected = chi2_contingency(contingency_table)
    p_values.append(p)
    
    # Calculate disease odds ratios
    homo_ref = contingency_table.loc[0, 1] / contingency_table.loc[0, 0] if 0 in contingency_table.index else np.nan
    het = contingency_table.loc[1, 1] / contingency_table.loc[1, 0] if 1 in contingency_table.index else np.nan
    homo_alt = contingency_table.loc[2, 1] / contingency_table.loc[2, 0] if 2 in contingency_table.index else np.nan

    odds_ratios_het.append(het / homo_ref if homo_ref != 0 else np.nan)
    odds_ratios_homo_alt.append(homo_alt / homo_ref if homo_ref != 0 else np.nan)


'''
# For PART A and B

# Convert p_values to a DataFrame
results = pd.DataFrame({'SNP': snp_ids, 'p_value': p_values})

# (b) Count how many SNPs have p-values below 0.05
significant_snps = results[results['p_value'] < 0.05]
print(f"Number of SNPs with p-value < 0.05: {len(significant_snps)}")

# Expected number of SNPs below 0.05 by chance
expected_by_chance = 0.05 * len(snp_ids)
print(f"Expected number of SNPs with p-value < 0.05 by chance: {expected_by_chance}")
'''

'''
# For part C

results = pd.DataFrame({
    'SNP': snp_ids,
    'uncorrected_p_value': p_values,
    'odds_ratio_het': odds_ratios_het,
    'odds_ratio_homo_alt': odds_ratios_homo_alt
})

# Apply Bonferroni correction
results['corrected_p_value'] = results['uncorrected_p_value'] * num_snps

# Filter significant SNPs after Bonferroni correction
significant_snps = results[results['corrected_p_value'] < 0.05]

# Display the results
print(significant_snps[['SNP', 'uncorrected_p_value', 'corrected_p_value', 'odds_ratio_het', 'odds_ratio_homo_alt']])

'''





