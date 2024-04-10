import argparse
import os
import pandas as pd
import scipy.stats as stats
import numpy as np
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

########################################################
# Parse args
########################################################

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--species', type=str, choices=['hg38', 'hg37', 'mm9', 'mm10','TAIR10'], help='species of genome')
parser.add_argument('-u', '--up_genes', type=str, help='name of up-regulated genes file')
parser.add_argument('-d', '--down_genes', type=str, help='name of down-regulated genes file')
parser.add_argument('-p', '--peak_file', type=str, help='name of region-covered genes file')
parser.add_argument('-dst', '--distance', type=int, choices=[1000,10_000,25_000,50_000, 100_000,500_000, 1_000_000], default=1_000_000,
                    help='maximum regulation distance to trace; defaults to 1000000')
parser.add_argument('-udst', '--user_distance_file', type=str, help='name of user distance file')
parser.add_argument('-uid', '--uid', type=int, help='save file for different users')
parser.add_argument('-color', '--color', type=str, help='color of the figure', default = "EB8A7A1E3758")
parser.add_argument('-title', '--title', type=str, help='title of the result figure', default = "Region Associated DEGs")
parser.add_argument('-test', '--test', type=str, choices=['hypergeom','binomial'], help='test to be used in RAD',default="hypergeom")

args = parser.parse_args()

#########################################################

def asterisk(x:float):
    if 0.01 < x < 0.05:
        return "*"
    elif 0.001 < x < 0.01:
        return "**"
    elif 0.0001 < x < 0.001:
        return "***"
    elif x < 0.0001:
        return "****"
    else:
        return False

#########################################################
# Genome info
#########################################################

GENOME_PATHS = {
    'hg38': './refdata/Homo_sapiens.GRCh38.97_pcg_chr.txt',
    'hg37': './refdata/Homo_sapiens.GRCh37.75_pcg_chr.txt',
    'mm9': './refdata/Mus_musculus.NCBIM37.67_pcg_chr.txt',
    'mm10': './refdata/Mus_musculus.GRCm38.97_pcg_chr.txt',
    'TAIR10': './refdata/Arabidopsis_thaliana.TAIR10.47_pcg_chr.txt',
}
GENOME_SIZES = {
    'hg38': 19957,
    'hg37': 22810,
    'mm9': 26873,
    'mm10': 21823,
    'TAIR10':27185,
}
genome_path = GENOME_PATHS[args.species]
genome_size = GENOME_SIZES[args.species]

##########################################################

##########################################################
# Read up/dw input
##########################################################
with open(args.up_genes) as up_f, open(args.down_genes) as dw_f:
    up_lines = up_f.read().splitlines()
    dw_lines = dw_f.read().splitlines()

first_line = up_lines[0]
if first_line in ['SYMBOL', 'ENSEMBL']: # has header
    up = up_lines[1:]
    dw = dw_lines[1:]
else:
    up = up_lines
    dw = dw_lines
up_origin_count = len(up)
dw_origin_count = len(dw)
##########################################################

##########################################################
# Identify which of ENSEMBL/SYMBOL is used
##########################################################
human_p = re.compile(r'^ENSG\d+')
mouse_p = re.compile(r'^ENSMUSG\d+')
tair_p = re.compile(r'^AT[0-9|M|C]G\d+')

is_ensembl = (first_line == 'ENSEMBL') or tair_p.match(first_line) or human_p.match(first_line) or mouse_p.match(first_line)
###########################################################

###########################################################
# Filter protein encoding genes
###########################################################
ref_pcg = pd.read_csv(genome_path,sep='\t',header = None)
if is_ensembl:
    ref_pcg = frozenset(ref_pcg.iloc[:,6])
else:
    ref_pcg = frozenset(ref_pcg.iloc[:,3])
is_pcg = lambda x: x in ref_pcg
up_pcg = filter(is_pcg, up) # iter, not list
dw_pcg = filter(is_pcg, dw) # iter, not list

###########################################################

color1 = args.color[0:6]
color2 = args.color[6:12]
title = args.title

###########################################################
# Generate peak-covered genes
###########################################################
with open(args.peak_file) as f:
    peak_first_element = f.readline().split("\t")[0]
if peak_first_element.startswith('chr'):
    txt_format = """awk '{{ print $1"\\t"(int(($2+$3)/2))"\\t"(int(($2+$3)/2)+1) }}' {0} \\
                | awk '{{ print $1"\\t"($2-{1}-1000)"\\t"($2+{1}+1000)"\\t"$3 }}' - \\
                | awk '{{ if ($2<0)print $1"\\t1\\t"$3"\\t"$4;else print $0 }}' - \\
                | bedtools intersect -a {2} -b - -wb > """.format(args.peak_file, args.distance, genome_path)
else:
    txt_format = """awk '{{ print "chr"$1"\\t"(int(($2+$3)/2))"\\t"(int(($2+$3)/2)+1) }}' {0} \\
                | awk '{{ print $1"\\t"($2-{1}-1000)"\\t"($2+{1}+1000)"\\t"$3 }}' - \\
                | awk '{{ if ($2<0)print $1"\\t1\\t"$3"\\t"$4;else print $0 }}' - \\
                | bedtools intersect -a {2} -b - -wb > """.format(args.peak_file, args.distance, genome_path)

###########################################################

covered_gene_path = "./static/file/{}-{}.txt".format('RAD_genename_distance', args.uid)
txt_format += covered_gene_path
# print("txt_format",txt_format)
os.system(txt_format)
cov = pd.read_table(covered_gene_path, header=None)
# print(cov.head())
# rebuild cov
cov_plus = cov.loc[cov.iloc[:, 5] == '+']
cov_minus = cov.loc[cov.iloc[:, 5] == '-']
# print('cov_plus')
# print(cov_plus.head())
# print("cov_minus")
# print(cov_minus.head())

cov_TSSdis = (cov_plus.iloc[:, 12] - cov_plus.iloc[:, 1]).append(cov_minus.iloc[:, 12] - cov_minus.iloc[:, 2])
cov = pd.DataFrame({
    'ENSEMBL': cov_plus.iloc[:, 6].append(cov_minus.iloc[:, 6]),
    'SYMBOL': cov_plus.iloc[:, 3].append(cov_minus.iloc[:, 3]),
})


# distance

LOWER_BOUNDS = {
    1_000: [-1000,-500,-200,-100,-50,0,
            50,100,200,500],
    10_000: [-10_000,-5_000,-2_000,-1_000,-500,-100,0,
            100,500,1000,2_000,5_000],
    25_000: [-25_000,-10_000,-5_000,-2_000,-1_000,-500,-100,0,
            100,500,1000,2_000,5_000,10_000],
    50_000: [-50_000,-25_000,-10_000,-5_000,-2_000,-1_000,-500,-100,0,
            100,500,1000,2_000,5_000,10_000,25_000],
    100_000: [-100_000,-50_000,-25_000,-10_000,-5_000,-2_000,-1_000,-500,-100,0,
            100,500,1000,2_000,5_000,10_000,25_000,50_000],
    500_000: [-500_000, -200_000, -100_000, -50_000, -25_000, -10_000, -5_000, -2_000,-1_000,0,
            1_000,2_000,5_000, 10_000, 25_000,50_000, 100_000, 200_000],
    1_000_000: [-1000_000, -500_000, -200_000, -100_000, -50_000, -25_000, -10_000, -5_000,-2_000, -1_000, 0,
            1_000,2_000,5_000, 10_000, 25_000,50_000, 100_000, 200_000, 500_000],
}
UPPER_BOUNDS = {
    1_000: [-500,-200,-100,-50,0,
            50,100,200,500,1000],
    10_000: [-5_000,-2_000,-1_000,-500,-100,0,
            100,500,1000,2_000,5_000,10_000],
    25_000: [-10_000,-5_000,-2_000,-1_000,-500,-100,0,
            100,500,1000,2_000,5_000,10_000,25_000],
    50_000: [-25_000,-10_000,-5_000,-2_000,-1_000,-500,-100,0,
            100,500,1000,2_000,5_000,10_000,25_000,50_000],
    100_000: [-50_000,-25_000,-10_000,-5_000,-2_000,-1_000,-500,-100,0,
            100,500,1000,2_000,5_000,10_000,25_000,50_000,100_000],
    500_000: [-200_000, -100_000, -50_000, -25_000, -10_000, -5_000, -2_000,-1_000,0,
            1_000,2_000,5_000, 10_000, 25_000,50_000, 100_000, 200_000, 500_000],
    1_000_000: [-500_000, -200_000, -100_000, -50_000, -25_000, -10_000, -5_000,-2_000, -1_000, 0,
            1_000,2_000,5_000, 10_000, 25_000,50_000, 100_000, 200_000, 500_000, 1000_000],
}
if args.user_distance_file:
    with open(args.user_distance_file) as f:
        temp_user_distance = f.read()
    temp_array = temp_user_distance.replace('#','').split(',')
    user_distance = [int(float(x) * 1000) for x in temp_array]
    user_distance_r = user_distance[::-1]
    neg_user_distance = [-1 * x for x in user_distance][:-1]
    user_dst = neg_user_distance + user_distance_r
    lower_bounds = user_dst[:-1]
    upper_bounds = user_dst[1:]
else:
    lower_bounds = LOWER_BOUNDS[args.distance]
    upper_bounds = UPPER_BOUNDS[args.distance]

# print(lower_bounds)


total_counts, up_counts, dw_counts, up_genes, dw_genes, up_dists, dw_dists = [], [], [], [], [], [], []
up_tbl = frozenset(up_pcg)
dw_tbl = frozenset(dw_pcg)
up_pcg_len = len(up_tbl) # duplicates removed
dw_pcg_len = len(dw_tbl)
# up_pcg = list(up_pcg)
# dw_pcg = list(dw_pcg)
# up_pcg_len = len(up_pcg)
# dw_pcg_len = len(dw_pcg)

assert (len(lower_bounds) == len(upper_bounds))
for low, high in zip(lower_bounds, upper_bounds):
    new_record = cov.loc[(cov_TSSdis < high) & (cov_TSSdis >= low)]
    total_counts.append(len(new_record))
    if is_ensembl:
        is_up = new_record['ENSEMBL'].apply(lambda x: x in up_tbl)
        is_dw = new_record['ENSEMBL'].apply(lambda x: x in dw_tbl)
    else:
        is_up = new_record['SYMBOL'].apply(lambda x: x in up_tbl)
        is_dw = new_record['SYMBOL'].apply(lambda x: x in dw_tbl)

    new_up_rec = new_record.loc[is_up]
    new_dw_rec = new_record.loc[is_dw]
    up_counts.append(len(new_up_rec))
    dw_counts.append(len(new_dw_rec))
    up_genes.append(new_up_rec)
    dw_genes.append(new_dw_rec)

    dist = "{:g}kb".format(low / 1000)
    up_dists.extend([dist, ] * up_counts[-1])
    dw_dists.extend([dist, ] * dw_counts[-1])

up_degs = pd.concat(up_genes).reset_index(drop=True)
up_degs['type'] = ['up', ] * len(up_dists)
up_degs['dist'] = up_dists
dw_degs = pd.concat(dw_genes).reset_index(drop=True)
dw_degs['type'] = ['dw', ] * len(dw_dists)
dw_degs['dist'] = dw_dists
degs = up_degs.append(dw_degs)

genenames_path = "./static/file/{}-{}.txt".format('RAD_genename_distance', args.uid)
degs.to_csv(genenames_path, sep='\t', index=False)

control = {
    'up': up_pcg_len,
    'dw': dw_pcg_len,
    'total': genome_size,
    'type': 'control',
    'dist': 'control',
}

dist_bins = ["{:g}kb to {:g}kb".format(low / 1000, high / 1000) for low, high in zip(lower_bounds, upper_bounds)]


if args.test == "hypergeom":
    counts = pd.DataFrame({
        'up': up_counts,
        'dw': dw_counts,
        'total': total_counts,
        'type': 'covered_genes',
        'dist': dist_bins,
        'pvalue_up': [stats.hypergeom.sf(up_count, control['total'], control['up'], total_count)
                      for up_count, total_count in zip(up_counts, total_counts)],
        'pvalue_dw': [stats.hypergeom.sf(dw_count, control['total'], control['dw'], total_count)
                      for dw_count, total_count in zip(dw_counts, total_counts)],
    })

    
elif args.test == "binomial":
    p1=control['up']/control['total']
    p2=control['dw']/control['total']

    counts = pd.DataFrame({
        'up': up_counts,
        'dw': dw_counts,
        'total': total_counts,
        'type': 'covered_genes',
        'dist': dist_bins,
        'pvalue_up': [stats.binom.sf(up_count, total_count, p1) for up_count, total_count in zip(up_counts, total_counts)],
        'pvalue_dw': [stats.binom.sf(dw_count, total_count, p2) for dw_count, total_count in zip(dw_counts, total_counts)]
    })
else:
    raise TypeError("Expected either hypergeometric test or binomial test")


pvalue_path = "./static/file/{}-{}.txt".format('RAD_genecount_pvalue', args.uid)
counts.to_csv(pvalue_path, sep='\t', index=False)

plot_up = (counts['up'] / counts['total']) / (control['up'] / control['total'])
plot_dw = (counts['dw'] / counts['total']) / (control['dw'] / control['total'])
plot_title = title

up_fail = (up_origin_count - up_pcg_len) / up_origin_count
dw_fail = (dw_origin_count - dw_pcg_len) / dw_origin_count
with open("./static/file/{}-{}.txt".format("test", args.uid),'w') as f:
    f.write("{}/{}, {:.2%}\t{}/{}, {:.2%}".format(
        up_origin_count - up_pcg_len, up_origin_count, up_fail,
        dw_origin_count - dw_pcg_len, dw_origin_count, dw_fail))

pdf = PdfPages("./static/file/{}-{}.pdf".format("Region Associated DEGs", args.uid))
plt.rcParams['font.sans-serif'] = 'Arial'
plt.figure(plot_title,figsize=(8,5),dpi = 300)
plt.title(plot_title)
plt.xlabel('Distance to TSS')
plt.ylabel('Observed/Expected')
plt.grid(linestyle=':', axis='y')
plt.tick_params(labelsize=10)
bottom, top = plt.ylim()

x = np.arange(len(lower_bounds))
up_bar = plt.bar(x - 0.2, plot_up, 0.4, label='up', color="#" + color1, align='center')
dw_bar = plt.bar(x + 0.2, plot_dw, 0.4, label='down', color="#" + color2, align='center')

for i, upp, dwp in zip(np.arange(len(lower_bounds)), counts['pvalue_up'], counts['pvalue_dw']):
    ast = asterisk(upp)
    if ast:
        plt.text(x=i-0.4,y=plot_up[i]+0.05,s=ast,size=8)
    ast = asterisk(dwp)
    if ast:
        plt.text(x=i,y=plot_dw[i]+0.05,s=ast,size=8)



plt.axhline(y=1, ls='--')
plt.xticks(x, dist_bins, rotation=-45,ha="left")
plt.legend()
plt.tight_layout()
# plt.show()

pdf.savefig()
plt.savefig("./static/file/{}-{}.png".format("Region Associated DEGs", args.uid))
plt.savefig("./static/file/{}-{}.svg".format("Region Associated DEGs", args.uid))
plt.close()
pdf.close()
