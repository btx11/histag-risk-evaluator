#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
動作確認・デモ用の de novo 酵素 + 基質 (SUB) + 補酵素 (NAD / FAD / SAM) ドッキング PDB ファイル生成スクリプト
"""

import math
from pathlib import Path
from typing import List, Tuple


def _format_atom_line(
    record: str,
    serial: int,
    atom_name: str,
    resname: str,
    chain: str,
    resseq: int,
    x: float,
    y: float,
    z: float,
    bfactor: float,
    element: str,
) -> str:
    name_field = f" {atom_name:<3s}" if len(atom_name) < 4 else f"{atom_name:<4s}"
    return (
        f"{record:<6s}{serial:5d} {name_field:4s} {resname:>3s} {chain:1s}"
        f"{resseq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{bfactor:6.2f}"
        f"          {element:>2s}\n"
    )


def build_synthetic_docked_enzyme(
    out_path: Path,
    n_term_pos: Tuple[float, float, float],
    c_term_pos: Tuple[float, float, float],
    bury_n_term: bool,
    bury_c_term: bool,
    substrate_center: Tuple[float, float, float],
    cofactor_center: Tuple[float, float, float],
    cofactor_name: str = "NAD",
    substrate_name: str = "SUB",
) -> None:
    """
    指定した N末・C末座標と露出状態、基質・補酵素ドッキング座標を持つテスト用 PDB ファイルを生成する。
    """
    lines: List[str] = [
        f"HEADER    DE NOVO ENZYME DOCKING COMPLEX          15-NOV-26   {out_path.stem[:4].upper()}\n",
        f"TITLE     SYNTHETIC DE NOVO ENZYME WITH DOCKED {substrate_name} AND {cofactor_name}\n",
    ]
    serial = 1
    sequence = [
        "MET", "GLY", "SER", "ALA", "GLU", "LEU", "LYS", "VAL", "ILE", "ALA",
        "ARG", "ASP", "TYR", "PHE", "TRP", "HIS", "ASN", "GLN", "PRO", "GLY",
        "ALA", "LEU", "GLU", "LYS", "VAL", "SER", "THR", "ILE", "ALA", "ARG",
        "ASP", "LEU", "TYR", "PHE", "ALA", "GLU", "LYS", "GLY", "SER", "ALA",
    ]
    n_res = len(sequence)

    # 1. 酵素鎖 (Chain A) の座標生成
    for i, resname in enumerate(sequence, start=1):
        t = (i - 1) / float(n_res - 1)
        if i <= 3:
            # N末端付近
            frac = (i - 1) / 3.0
            cx = n_term_pos[0] * (1.0 - frac) + 0.0 * frac
            cy = n_term_pos[1] * (1.0 - frac) + 0.0 * frac
            cz = n_term_pos[2] * (1.0 - frac) + 0.0 * frac
        elif i >= n_res - 2:
            # C末端付近
            frac = (n_res - i) / 3.0
            cx = c_term_pos[0] * (1.0 - frac) + 0.0 * frac
            cy = c_term_pos[1] * (1.0 - frac) + 0.0 * frac
            cz = c_term_pos[2] * (1.0 - frac) + 0.0 * frac
        else:
            # 酵素コア (らせんバンドル風配置)
            angle = t * 6.0 * math.pi
            radius = 9.5 + 2.5 * math.sin(t * 2.0 * math.pi)
            cx = radius * math.cos(angle)
            cy = radius * math.sin(angle)
            cz = (t - 0.5) * 22.0

        # もし N末 または C末 を埋没させる設定なら、周囲のコア残基を末端の周囲に密に配置する
        if bury_n_term and 10 <= i <= 22:
            ang = (i - 10) * (2.0 * math.pi / 12.0)
            r_shell = 3.6 if i % 2 == 0 else 4.8
            cx = n_term_pos[0] + r_shell * math.cos(ang)
            cy = n_term_pos[1] + r_shell * math.sin(ang)
            cz = n_term_pos[2] + ((i % 3) - 1) * 2.6
        elif bury_c_term and 15 <= i <= 27:
            ang = (i - 15) * (2.0 * math.pi / 12.0)
            r_shell = 3.6 if i % 2 == 0 else 4.8
            cx = c_term_pos[0] + r_shell * math.cos(ang)
            cy = c_term_pos[1] + r_shell * math.sin(ang)
            cz = c_term_pos[2] + ((i % 3) - 1) * 2.6

        bfac = 88.5 if 4 <= i <= n_res - 3 else 76.0
        atom_offsets = [
            ("N", -1.2, 0.4, 0.0, "N"),
            ("CA", 0.0, 0.0, 0.0, "C"),
            ("C", 1.3, 0.5, 0.0, "C"),
            ("O", 1.4, 1.7, 0.0, "O"),
        ]
        if resname != "GLY":
            atom_offsets.append(("CB", -0.3, -1.4, 0.6, "C"))
            atom_offsets.append(("CG", -0.6, -2.5, 1.2, "C"))

        for aname, dx, dy, dz, el in atom_offsets:
            lines.append(
                _format_atom_line(
                    "ATOM", serial, aname, resname, "A", i,
                    cx + dx, cy + dy, cz + dz, bfac, el
                )
            )
            serial += 1

    lines.append("TER\n")

    # 2. 補酵素 (Cofactor, Chain A, resseq 501)
    cof_atoms = [
        ("C1", 0.0, 0.0, 0.0, "C"), ("C2", 1.4, 0.2, 0.0, "C"),
        ("N1", 2.1, 1.3, 0.3, "N"), ("C3", 1.5, 2.5, 0.1, "C"),
        ("O1", 2.2, 3.5, 0.2, "O"), ("P1", -1.5, -0.7, 0.4, "P"),
        ("O2", -2.5, 0.3, 0.7, "O"), ("C4", -1.8, -1.9, -0.8, "C"),
    ]
    for aname, dx, dy, dz, el in cof_atoms:
        lines.append(
            _format_atom_line(
                "HETATM", serial, aname, cofactor_name, "A", 501,
                cofactor_center[0] + dx, cofactor_center[1] + dy, cofactor_center[2] + dz,
                82.0, el,
            )
        )
        serial += 1

    # 3. 基質 (Substrate, Chain A, resseq 502)
    sub_atoms = [
        ("C1", 0.0, 0.0, 0.0, "C"), ("C2", 1.4, 0.0, 0.2, "C"),
        ("C3", 2.1, 1.2, 0.0, "C"), ("O1", 3.3, 1.2, 0.1, "O"),
        ("N1", 1.4, 2.4, -0.3, "N"), ("C4", 0.0, 2.3, -0.2, "C"),
    ]
    for aname, dx, dy, dz, el in sub_atoms:
        lines.append(
            _format_atom_line(
                "HETATM", serial, aname, substrate_name, "A", 502,
                substrate_center[0] + dx, substrate_center[1] + dy, substrate_center[2] + dz,
                79.0, el,
            )
        )
        serial += 1

    lines.append("END\n")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("".join(lines), encoding="utf-8")


def generate_all_samples(sample_dir: Path = Path("./sample_pdbs")) -> List[Path]:
    sample_dir.mkdir(parents=True, exist_ok=True)
    p1 = sample_dir / "denovo_enz01_N_buried_near_active.pdb"
    build_synthetic_docked_enzyme(
        out_path=p1,
        # N末が基質ポケットのすぐ近く (約 7.5 Å) かつ埋没、C末は遠く (約 28 Å) かつ溶媒露出
        n_term_pos=(4.0, 3.0, 2.0),
        c_term_pos=(-18.0, -16.0, 14.0),
        bury_n_term=True,
        bury_c_term=False,
        substrate_center=(8.0, 6.0, 4.0),
        cofactor_center=(5.0, 8.0, 6.0),
        cofactor_name="NAD",
        substrate_name="SUB",
    )

    p2 = sample_dir / "denovo_enz02_C_near_active.pdb"
    build_synthetic_docked_enzyme(
        out_path=p2,
        # N末は遠く (約 27 Å) かつ溶媒露出、C末が基質・補酵素のすぐ近く (約 8.0 Å)
        n_term_pos=(-17.0, 16.0, -13.0),
        c_term_pos=(6.0, -4.0, 3.0),
        bury_n_term=False,
        bury_c_term=False,
        substrate_center=(9.0, -7.0, 5.0),
        cofactor_center=(7.0, -9.0, 2.0),
        cofactor_name="FAD",
        substrate_name="SUB",
    )

    p3 = sample_dir / "denovo_enz03_both_exposed_distal.pdb"
    build_synthetic_docked_enzyme(
        out_path=p3,
        # N末・C末ともに活性部位から十分離れており (> 24 Å)、かつ両方とも溶媒露出
        n_term_pos=(-16.0, -15.0, -12.0),
        c_term_pos=(18.0, 16.0, 14.0),
        bury_n_term=False,
        bury_c_term=False,
        substrate_center=(2.0, -2.0, 1.0),
        cofactor_center=(-2.0, 2.0, -1.0),
        cofactor_name="SAM",
        substrate_name="SUB",
    )
    return [p1, p2, p3]


if __name__ == "__main__":
    paths = generate_all_samples()
    for p in paths:
        print(f"Generated sample PDB: {p}")
