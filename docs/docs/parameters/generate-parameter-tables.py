import xml.etree.ElementTree as ET


groups = {
    'w90': [
        'disentanglement',
        'job',
        'plot',
        'system',
        'transport',
        'wannierise',
    ],
    'postw90': [
        'berry',
        'boltzwann',
        'dos',
        'geninterp',
        'global',
        'gyrotropic',
        'kpath',
        'kslice',
    ],
}

tree = ET.parse('parameters.xml')
root = tree.getroot()

for tool in ['w90', 'postw90']:
    for group in groups[tool]:
        parameters = root.findall(f'./parameter[@tool="{tool}"][@group="{group}"]')
        with open(f'{tool}-{group}-parameters.csv', 'w') as fp:
            print('Keyword,Type,Description', file=fp)
            for parameter in parameters:
                name = parameter.find('name')
                type = parameter.find('type')
                description = parameter.find('description')
                if 'optional_prefix' in name.attrib:
                    print(f'[{name.attrib["optional_prefix"]}]{name.text},', end='', file=fp)
                else:
                    print(f'{name.text},', end='', file=fp)
                print(f'{type.text},', end='', file=fp)
                print(f'"{description.text}"', file=fp)
