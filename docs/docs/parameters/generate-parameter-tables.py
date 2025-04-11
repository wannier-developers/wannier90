import xml.etree.ElementTree as ET


groups = {
    "w90": [
        "disentanglement",
        "job",
        "plot",
        "system",
        "transport",
        "wannierise",
    ],
    "postw90": [
        "berry",
        "boltzwann",
        "dos",
        "geninterp",
        "global",
        "gyrotropic",
        "kpath",
        "kslice",
    ],
}

tree = ET.parse("parameters.xml")
root = tree.getroot()

for tool in ["w90", "postw90"]:
    tool_parameters = root.findall(f'./parameter[@tool="{tool}"]')
    for group in groups[tool]:
        parameters = [p for p in tool_parameters if any([g.text == group for g in p.find("groups").findall("group")])]
        with open(f"{tool}-{group}-parameters.csv", "w") as fp:
            print("Keyword,Type,Description", file=fp)
            for parameter in parameters:
                name = parameter.find("name")
                types = parameter.find("types")
                units = None
                if types:
                    type_text = '"' + ",".join([t.text for t  in types.findall("type")]) + '"'
                else:
                    type = parameter.find("type")
                    units = type.get("units")
                    type_text = type.text
                description = parameter.find("description")
                if "optional_prefix" in name.attrib:
                    print(f'[{name.attrib["optional_prefix"]}]{name.text},', end="", file=fp)
                else:
                    print(f"{name.text},", end="", file=fp)
                print(f"{type_text},", end="", file=fp)
                description_text = description.text
                if units and "$" not in units:
                    description_text += f" ({units})"
                print(f'"{description_text}"', file=fp)
