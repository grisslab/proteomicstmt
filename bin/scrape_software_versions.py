#!/usr/bin/env python

from __future__ import print_function
from collections import OrderedDict
import re

openms_version_regex = r"([0-9][.][0-9][.][0-9])"

regexes = {
    'nf-core/proteomicsTMT': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'ThermorawfileParser': ['v_thermorawfileparser.txt', r"(\S+)"],
	'IsobaricAnalyzer': ['v_isobaricanalyzer.txt', openms_version_regex],
    'FileConverter': ['v_fileconverter.txt', openms_version_regex],
    'DecoyDatabase': ['v_decoydatabase.txt', openms_version_regex],
    'MSGFPlusAdapter': ['v_msgfplusadapter.txt', openms_version_regex],
    'MSGFPlus': ['v_msgfplus.txt', r"\(([^v)]+)\)"],
    'CometAdapter': ['v_cometadapter.txt', openms_version_regex],
    'Comet': ['v_comet.txt', r"\"(.*)\""],
    'PeptideIndexer': ['v_peptideindexer.txt', openms_version_regex],
    'PSMFeatureExtractor': ['v_psmfeatureextractor.txt', openms_version_regex],
    'PercolatorAdapter': ['v_percolatoradapter.txt', openms_version_regex],
    'Percolator': ['v_percolator.txt', r"([0-9].[0-9]{2}.[0-9])"],
    'IDFilter': ['v_idfilter.txt', openms_version_regex],
    'IDScoreSwitcher': ['v_idscoreswitcher.txt', openms_version_regex],
    'FalseDiscoveryRate': ['v_falsediscoveryrate.txt', openms_version_regex],
    'IDPosteriorErrorProbability': ['v_idposteriorerrorprobability.txt', openms_version_regex],
    'IDMapper': ['v_idmapper.txt', openms_version_regex],
	'FileMerger': ['v_filemerger.txt', openms_version_regex],
	'Epifany': ['v_epifany.txt', openms_version_regex],
	'IDConflictResolver': ['v_idconflictresolver', openms_version_regex],
	'ProteinQuantifier': ['v_proteinquantifier.txt', openms_version_regex],
	'MSstatsConverter': ['v_msstatsconverter.txt', openms_version_regex],
    'MSstatsTMT': ['v_msstats_tmt.txt', r"(\S+)"],
	'pmultiqc': ['v_pmultiqc.txt', openms_version_regex]
}
results = OrderedDict()
results['nf-core/proteomicsTMT'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['ThermorawfileParser'] = '<span style="color:#999999;\">N/A</span>'
results['IsobaricAnalyzer'] = '<span style="color:#999999;\">N/A</span>'
results['FileConverter'] = '<span style="color:#999999;\">N/A</span>'
results['DecoyDatabase'] = '<span style="color:#999999;\">N/A</span>'
results['MSGFPlusAdapter'] = '<span style="color:#999999;\">N/A</span>'
results['MSGFPlus'] = '<span style="color:#999999;\">N/A</span>'
results['CometAdapter'] = '<span style="color:#999999;\">N/A</span>'
results['Comet'] = '<span style="color:#999999;\">N/A</span>'
results['PeptideIndexer'] = '<span style="color:#999999;\">N/A</span>'
results['PSMFeatureExtractor'] = '<span style="color:#999999;\">N/A</span>'
results['PercolatorAdapter'] = '<span style="color:#999999;\">N/A</span>'
results['Percolator'] = '<span style="color:#999999;\">N/A</span>'
results['IDFilter'] = '<span style="color:#999999;\">N/A</span>'
results['IDScoreSwitcher'] = '<span style="color:#999999;\">N/A</span>'
results['FalseDiscoveryRate'] = '<span style="color:#999999;\">N/A</span>'
results['IDPosteriorErrorProbability'] = '<span style="color:#999999;\">N/A</span>'
results['IDMapper'] = '<span style="color:#999999;\">N/A</span>'
results['FileMerger'] = '<span style="color:#999999;\">N/A</span>'
results['Epifany'] = '<span style="color:#999999;\">N/A</span>'
results['IDConflictResolver'] = '<span style="color:#999999;\">N/A</span>'
results['ProteinQuantifier'] = '<span style="color:#999999;\">N/A</span>'
results['MSstatsConverter'] = '<span style="color:#999999;\">N/A</span>'
results['MSstatsTMT'] = '<span style="color:#999999;\">N/A</span>'
results['pmultiqc'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    try:
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))
    except IOError:
        results[k] = False

# Remove software set to false in results
for k in list(results):
    if not results[k]:
        del results[k]

# Dump to YAML
print(
    """
id: 'software_versions'
section_name: 'nf-core/proteomicsTMT Software Versions'
section_href: 'https://github.com/nf-core/proteomicsTMT'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
    """
)
for k, v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k, v))
print("    </dl>")

# Write out regexes as csv file:
with open("software_versions.csv", "w") as f:
    for k, v in results.items():
        f.write("{}\t{}\n".format(k, v))
