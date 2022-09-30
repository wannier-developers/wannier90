#!/usr/bin/env python3
"""
A simple script to extract pseudo-atomic orbitals in UPF format
to a dat file.

The dat file can be read by ``pw2wannier90.x`` to compute
Amn matrices using these external projectors.
"""
from pathlib import Path
from lxml import etree
import numpy as np
import matplotlib.pyplot as plt
import sys

from collections import namedtuple

PSWFC = namedtuple("PSWFC", ["label", "data"])


def element2array(xml_element):
    txt = xml_element.text
    txt = txt.replace("\n", " ")
    arr = np.fromstring(txt, sep=" ")
    return arr


def read_upf_rmesh(filename):
    tree = etree.parse(filename)
    root = tree.getroot()
    rmesh = element2array(root.find("./PP_MESH/PP_R"))
    # print(rmesh)
    return rmesh


def read_upf_xpath(filename, xpath):
    if not isinstance(xpath, str):
        raise ValueError(f"only accept str type for parameter xpath: `{xpath}`")

    tree = etree.parse(filename)
    root = tree.getroot()
    res = {}
    for element in root.xpath(xpath):
        name = element.tag
        print(name)
        assert name != ""
        assert name not in res
        res[name] = element2array(element)
    return res


def tag2xpath(tags):
    if isinstance(tags, str):
        tags = [tags]

    if not isinstance(tags, (list, tuple)):
        raise Exception(f"wrong argument type: {tags}")

    tag2xpath = {
        "PP_NLCC": "PP_NLCC",
        "PP_LOCAL": "PP_LOCAL",
        "PP_NONLOCAL/PP_BETA": "PP_NONLOCAL/*[starts-with(local-name(), 'PP_BETA.')]",
        "PP_PSWFC/PP_CHI": "PP_PSWFC/*[starts-with(local-name(), 'PP_CHI.')]",
        "PP_FULL_WFC/PP_AEWFC": "PP_FULL_WFC/*[starts-with(local-name(), 'PP_AEWFC.')]",
        "PP_FULL_WFC/PP_PSWFC": "PP_FULL_WFC/*[starts-with(local-name(), 'PP_PSWFC.')]",
        "PP_RHOATOM": "PP_RHOATOM",
        "PP_PAW/PP_AE_NLCC": "PP_PAW/PP_AE_NLCC",
        "PP_PAW/PP_AE_VLOC": "PP_PAW/PP_AE_VLOC",
    }

    xpaths = []
    for tag in tags:
        xpaths.append(tag2xpath[tag])
    return xpaths


def read_vanderbilt_rmesh(filename):
    """Read UPF in format like 'ti_pbe_v1.4.uspp.F.UPF'."""
    with open(filename) as fil:
        lines = fil.readlines()
        idx_start = lines.index("  <PP_R>\n")
        idx_stop = lines.index("  </PP_R>\n")
        # print(idx_start, idx_stop)
        # print(lines[idx_start+1:idx_stop])
        lines = " ".join(lines[idx_start + 1 : idx_stop])
        lines = lines.replace("\n", " ")
        rmesh = np.fromstring(lines, sep=" ")
        return rmesh


def read_vanderbilt_pswfc(filename):
    """Read UPF in format like 'ti_pbe_v1.4.uspp.F.UPF'."""
    with open(filename) as fil:
        lines = fil.readlines()
        idx_start = lines.index("<PP_PSWFC>\n")
        idx_stop = lines.index("</PP_PSWFC>\n")
        # print(idx_start, idx_stop)
        lines = lines[idx_start + 1 : idx_stop]
        # print(lines)
        pswfcs = []
        idx_last_wfc = 0
        for i, line in enumerate(lines):
            if "Wavefunction" in line or i == len(lines) - 1:
                idx_current_wfc = i
                if i == len(lines) - 1:
                    idx_current_wfc += 1
                if idx_current_wfc != idx_last_wfc:
                    label = lines[idx_last_wfc].split()[0]
                    data = lines[idx_last_wfc + 1 : idx_current_wfc]
                    data = " ".join(data)
                    data = data.replace("\n", " ")
                    data = np.fromstring(data, sep=" ")
                    pswfc = PSWFC(label, data)
                    idx_last_wfc = idx_current_wfc
                    pswfcs.append(pswfc)
        return pswfcs


def write_dat(upf_filename, dat_filename):
    rmesh = read_vanderbilt_rmesh(upf_filename)
    pswfcs = read_vanderbilt_pswfc(upf_filename)
    _write_dat(dat_filename, rmesh, pswfcs)


def _write_dat(dat_filename, rmesh, pswfcs):
    xmesh = np.log(rmesh)
    print(xmesh[1:] - xmesh[:-1])

    num_grid = len(rmesh)
    num_proj = len(pswfcs)

    map_l = {"s": 0, "p": 1, "d": 2, "f": 3}
    n_of_pswfc = lambda _: int(_.label[0])
    l_of_pswfc = lambda _: map_l[_.label[1].lower()]

    # First sort by l, then by n, in ascending order
    pswfcs = sorted(pswfcs, key=lambda _: (l_of_pswfc(_), n_of_pswfc(_)))
    # print(pswfcs)

    with open(dat_filename, "w") as fil:
        fil.write(f"{num_grid} {num_proj}\n")
        fil.write(" ".join([str(l_of_pswfc(_)) for _ in pswfcs]) + "\n")
        for i in range(num_grid):
            fil.write(f"{xmesh[i]:.15f}")
            fil.write(f" {rmesh[i]:.15f}")
            for j in range(num_proj):
                fil.write(f" {pswfcs[j].data[i]}")
            fil.write("\n")


def plot_vanderbilt_upf(filename):
    """For UPF like 'ti_pbe_v1.4.uspp.F.UPF'."""

    rmesh = read_vanderbilt_rmesh(filename)
    rmax = rmesh[0]
    pswfcs = read_vanderbilt_pswfc(filename)
    print(len(rmesh))
    for name, array in pswfcs:
        print(name, len(array))
        # exclude_indexes = [1,2,3,4]
        exclude_indexes = []
        if np.any([str(ex_ind) in name for ex_ind in exclude_indexes]):
            continue
        # limit xrange in plot
        ind = np.argwhere(np.abs(array) > 1e-3).flatten()[-1]
        # print(ind, rmesh[ind], arr[ind-5:ind+5])
        rmax = max(rmax, rmesh[min(ind + 3, len(rmesh) - 1)])
        plt.plot(rmesh, array, label=name)
    plt.axhline(y=0, linewidth=0.4, color="k", linestyle="--")
    plt.xlim(rmesh[0], rmax)
    plt.xlabel("r (Bohr)")
    plt.ylabel("PSWFC")
    plt.title(filename)
    plt.legend()
    plt.show()


def plot_upf(filename, tags):
    """plot components of UPF

    :param filename: [description]
    :type filename: [type]
    :param tags: supported tags, see the keys of tag2xpath
    :type tags: str or list or tuple
    """
    # is_pswfc = True
    is_pswfc = False

    xpaths = tag2xpath(tags)
    rmesh = read_upf_rmesh(filename)
    rmax = rmesh[0]
    for xp in xpaths:
        for name, array in read_upf_xpath(filename, xp).items():
            # exclude_indexes = [1,2,3,4]
            exclude_indexes = []
            if np.any([str(ex_ind) in name for ex_ind in exclude_indexes]):
                continue
            # limit xrange in plot
            ind = np.argwhere(np.abs(array) > 1e-3).flatten()[-1]
            # print(ind, rmesh[ind], arr[ind-5:ind+5])
            rmax = max(rmax, rmesh[min(ind + 3, len(rmesh) - 1)])
            if is_pswfc:
                name = "3S" if "1" in name else "3P"
            plt.plot(rmesh, array, label=name)
    plt.axhline(y=0, linewidth=0.4, color="k", linestyle="--")
    plt.xlim(rmesh[0], rmax)
    plt.xlabel("r (Bohr)")
    if is_pswfc:
        plt.ylabel("PSWFC")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: write_pdwf_projectors.py <upf_file>")
        sys.exit(1)

    filename = sys.argv[1]

    # Example usage:
    #
    # filename = 'Si.pbe-nl-kjpaw_psl.1.0.0.UPF'
    # filename = 'ti_pbe_v1.4.uspp.F.UPF'
    # plot_vanderbilt_upf(filename)

    # filename = 'Ga.pbe-dn-kjpaw_psl.1.0.0.UPF'
    # filename = "P.pbe-n-rrkjus_psl.1.0.0.UPF"
    #
    # plot_upf(filename, "PP_PSWFC/PP_CHI")
    # plot_upf(filename, ['PP_NLCC', 'PP_PAW/PP_AE_NLCC'])
    # plot_upf(filename, 'PP_LOCAL')
    # plot_upf(filename, 'PP_NONLOCAL/PP_BETA')
    plot_upf(filename, "PP_PSWFC/PP_CHI")
    # plot_upf(filename, 'PP_FULL_WFC/PP_AEWFC')
    # plot_upf(filename, 'PP_FULL_WFC/PP_PSWFC')
    # plot_upf(filename, ['PP_PSWFC/PP_CHI', 'PP_FULL_WFC/PP_PSWFC'])
    # plot_upf(filename, ['PP_FULL_WFC/PP_PSWFC', 'PP_FULL_WFC/PP_AEWFC'])
    # plot_upf(filename, 'PP_RHOATOM')
    # plot_upf(filename, 'PP_PAW/PP_AE_NLCC')
    # plot_upf(filename, ['PP_LOCAL', 'PP_PAW/PP_AE_VLOC'])

    write_dat(filename, f'{filename}.dat')
