import re
import xml.etree.ElementTree as ET

import pandas as pd


csv_files = [
    'w90-system-parameters.csv',
    'w90-job-parameters.csv',
    'w90-disentanglement-parameters.csv',
    'w90-wannierise-parameters.csv',
    'w90-plot-parameters.csv',
    'w90-transport-parameters.csv',
    'postw90-global-parameters.csv',
    'postw90-berry-parameters.csv',
    'postw90-dos-parameters.csv',
    'postw90-kpath-parameters.csv',
    'postw90-kslice-parameters.csv',
    'postw90-gyrotropic-parameters.csv',
    'postw90-boltzwann-parameters.csv',
    'postw90-geninterp-parameters.csv',
]

root = ET.Element('root')
tree = ET.ElementTree(root)

pattern = re.compile(r'(?P<tool>.+)-(?P<group>.+)-parameters.csv')

for csv_file in csv_files:
    match1 = pattern.match(csv_file)

    table = pd.read_csv(
            csv_file,
                header=0,
                names=('keyword', 'type', 'description'),
                converters={'keyword': str.strip, 'type': str.strip, 'description': str.strip}
            )

    for row in table.itertuples():
        parameter = ET.SubElement(root, 'parameter')
        parameter.set('tool', match1['tool'])
        parameter.set('group', match1['group'])
        # tool = ET.SubElement(parameter, 'tool')
        # group = ET.SubElement(parameter, 'group')
        name = ET.SubElement(parameter, 'name')
        # optional_prefix = ET.SubElement(parameter, 'optional_prefix')
        type = ET.SubElement(parameter, 'type')
        description = ET.SubElement(parameter, 'description')

        if (match2 := re.compile(r'^\[(?P<prefix>.+_)\]').match(row.keyword)) is not None:
            name.text = row.keyword.replace(f'[{match2["prefix"]}]', '')
            # optional_prefix.text = match2['prefix']
            name.set('optional_prefix', match2['prefix'])
        else:
            name.text = row.keyword
        # tool.text = match1['tool']
        # group.text = match1['group']
        type.text = row.type
        description.text = row.description

ET.indent(tree)
with open('parameters.xml', 'wb') as fp:
    tree.write(fp)
