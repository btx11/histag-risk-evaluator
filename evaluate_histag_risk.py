#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
de novo 酵素の N末端・C末端 Hisタグ付加リスク評価 & PyMOL 3方向バッチ作図スクリプト
================================================================================
ディレクトリ内の複数の酵素-リガンド（基質・補酵素）ドッキング立体構造 PDB ファイルを一括解析し、
1. PyMOL を用いて以下の配色で 3 方向（正面 0°、90°回転×2種）の立体構造図をバッチ作成:
   - 酵素（タンパク質鎖）: 青色 (blue)
   - リガンド（基質・補酵素の2種など）: 黄色 (yellow)
   - N末端: 赤色 (red)
   - C末端: 緑色 (green)
   ※ Windows PC 上のデスクトップ版 PyMOL および WSL2 (Ubuntu) 環境の両方に完全対応。
2. 構造生物学・タンパク質工学の先行研究（論文）に基づき、N末端・C末端それぞれの
   - リガンド（基質・補酵素）までの立体距離（最小重原子間距離 d_min および Cα-重心間距離）
   - 溶媒露出表面積 (SASA; Shrake & Rupley, 1973) と相対溶媒露出度 (RSA; Tien et al., 2013)
   - 末端主鎖ベクトルの活性部位に対する配向角 (Orientation angle)
   を算出し、Hisタグ付加リスクスコア (0–100) と推奨タグ位置を自動判定。
3. 算出された数値・判定結果・引用文献情報・3方向の構造画像を Excel (.xlsx) に自動集約。

主な引用文献 (References):
--------------------------
[1] Shrake, A., & Rupley, J. A. (1973). Environment and exposure to solvent of protein atoms.
    Lysozyme and insulin. Journal of Molecular Biology, 79(2), 351-371.
    DOI: 10.1016/0022-2836(73)90011-9
[2] Tien, M. Z., Meyer, A. G., Sydykova, D. K., Spielman, S. J., & Wilke, C. O. (2013).
    Maximum allowed solvent accessibilities of residues in proteins. PLOS ONE, 8(11), e80635.
    DOI: 10.1371/journal.pone.0080635
[3] Thornton, J. M., & Sibanda, B. L. (1983). Amino and carboxy-terminal regions in globular
    proteins. Journal of Molecular Biology, 167(2), 443-460.
    DOI: 10.1016/S0022-2836(83)80344-1
[4] Miller, S., Janin, J., Lesk, A. M., & Chothia, C. (1987). Interior and surface of
    monomeric proteins. Journal of Molecular Biology, 196(3), 641-656.
    DOI: 10.1016/0022-2836(87)90038-6
[5] Rost, B., & Sander, C. (1994). Conservation and prediction of solvent accessibility in
    protein families. Proteins: Structure, Function, and Bioinformatics, 20(3), 216-226.
    DOI: 10.1002/prot.340200303
[6] Carson, M., Johnson, D. H., McDonald, H., Brouillette, C., & DeLucas, L. J. (2007).
    His-tag impact on structure. Acta Crystallographica Section D: Biological Crystallography,
    63(3), 295-301. DOI: 10.1107/S0907444906052024
[7] Freydank, A. C., Brandt, W., & Dräger, B. (2008). Protein structure modeling indicates
    hexahistidine-tag interference with enzyme activity. Proteins: Structure, Function, and
    Bioinformatics, 72(1), 173-183. DOI: 10.1002/prot.21905
[8] Yeon, Y. J., Park, H. J., & Yoo, Y. J. (2014). Effect of His-tag location on the
    catalytic activity of 3-hydroxybutyrate dehydrogenase. Biotechnology and Bioprocess
    Engineering, 19(5), 798-802. DOI: 10.1007/s12257-014-0089-2
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import platform
import shutil
import struct
import subprocess
import sys
import tempfile
import zlib
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple
from xml.sax.saxutils import escape as xml_escape

# =============================================================================
# 1. 定数・構造生物学パラメータ定義 (Literature Constants)
# =============================================================================

# 標準20アミノ酸および一般的なプロトン化・修飾表記 -> 標準3文字コードへのマッピング
PROTEIN_RESIDUE_MAP: Dict[str, str] = {
    "ALA": "ALA", "ARG": "ARG", "ASN": "ASN", "ASP": "ASP", "CYS": "CYS",
    "GLN": "GLN", "GLU": "GLU", "GLY": "GLY", "HIS": "HIS", "ILE": "ILE",
    "LEU": "LEU", "LYS": "LYS", "MET": "MET", "PHE": "PHE", "PRO": "PRO",
    "SER": "SER", "THR": "THR", "TRP": "TRP", "TYR": "TYR", "VAL": "VAL",
    # よく使われるプロトン化・セレノメチオニン等のバリアント
    "MSE": "MET", "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "CYX": "CYS", "CSO": "CYS",
    "ASH": "ASP", "GLH": "GLU", "LYN": "LYS", "SEP": "SER", "TPO": "THR",
    "PTR": "TYR",
}

# Tien et al. (2013) PLOS ONE 8(11): e80635, Table 1
# Theoretical Maximum Allowed Solvent Accessibility (MaxASA, Å^2) for Gly-X-Gly tripeptides
TIEN_2013_MAX_ASA_THEORETICAL: Dict[str, float] = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLU": 223.0, "GLN": 225.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}

# Tien et al. (2013) Empirical MaxASA (Å^2)
TIEN_2013_MAX_ASA_EMPIRICAL: Dict[str, float] = {
    "ALA": 121.0, "ARG": 265.0, "ASN": 187.0, "ASP": 187.0, "CYS": 148.0,
    "GLU": 214.0, "GLN": 214.0, "GLY":  97.0, "HIS": 216.0, "ILE": 195.0,
    "LEU": 191.0, "LYS": 230.0, "MET": 203.0, "PHE": 228.0, "PRO": 154.0,
    "SER": 143.0, "THR": 163.0, "TRP": 264.0, "TYR": 255.0, "VAL": 165.0,
}

# Bondi (1964) / Chothia (1976) 原子ファンデルワールス半径 (Å)
VDW_RADII: Dict[str, float] = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80,
    "P": 1.80, "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98,
    "FE": 1.40, "ZN": 1.39, "MG": 1.73, "CA": 1.97, "MN": 1.40,
    "CU": 1.40, "CO": 1.40, "NI": 1.63, "SE": 1.90,
}
DEFAULT_VDW_RADIUS: float = 1.70
WATER_PROBE_RADIUS: float = 1.40  # Shrake & Rupley (1973) 水分子プローブ半径 1.4 Å

# 水分子および一般的な結晶化緩衝液・塩イオン（リガンドから除外するデフォルト一覧）
WATER_RESNAMES = {"HOH", "WAT", "TIP", "TIP3", "TP3", "SOL", "DOD", "D2O"}
BUFFER_AND_SALT_RESNAMES = {
    "NA", "CL", "K", "BR", "IOD", "LI", "RB", "CS",
    "SO4", "PO4", "GOL", "EDO", "PEG", "PGE", "DMS", "ACT", "ACY",
    "TRS", "EPE", "MES", "CIT", "MPD", "FMT", "NO3", "SCN", "URE", "bme", "BME",
}

# 代表的な補酵素 (Cofactor) の残基名リスト（基質と補酵素の自動判別用）
KNOWN_COFACTOR_RESNAMES = {
    # ニコチンアミド系
    "NAD", "NAI", "NAP", "NDP", "NAJ", "NHD", "NMN",
    # フラビン系
    "FAD", "FMN", "FNR", "FDA",
    # ヘム・ポルフィリン・鉄硫黄クラスター系
    "HEM", "HEC", "HEA", "HEB", "HEO", "FES", "SF4", "F3S", "F4S",
    # アデノシル・ヌクレオチド・エネルギー補酵素系
    "SAM", "SAH", "ATP", "ADP", "AMP", "ANP", "ACP", "GTP", "GDP", "CTP", "UTP",
    # ピリドキサール・チアミン・ビオチン・コエンザイムA・PQQ・コバラミン・プテリン系
    "PLP", "PMP", "TPP", "TDP", "THDP", "COA", "ACO", "DCA", "sca",
    "BIO", "BTN", "PQQ", "B12", "CNC", "COB", "MGD", "THF", "BH4", "H4B", "ASC",
    # 触媒金属補因子（オプションで補酵素扱いされる場合）
    "FE", "FE2", "ZN", "MN", "CU", "CO", "NI", "MG", "MO",
}


# =============================================================================
# 2. データクラス定義 (Data Structures)
# =============================================================================

@dataclass
class Atom:
    serial: int
    name: str
    altloc: str
    resname: str
    chain: str
    resseq: int
    icode: str
    x: float
    y: float
    z: float
    occupancy: float
    bfactor: float
    element: str
    record_type: str  # "ATOM" or "HETATM"
    sasa: float = 0.0

    @property
    def is_hydrogen(self) -> bool:
        if self.element.upper() in ("H", "D"):
            return True
        stripped = self.name.strip()
        return stripped.startswith("H") or (len(stripped) >= 2 and stripped[0].isdigit() and stripped[1] == "H")

    @property
    def vdw_radius(self) -> float:
        el = self.element.upper().strip()
        if el in VDW_RADII:
            return VDW_RADII[el]
        # 原子名から元素を推定
        name_clean = "".join([c for c in self.name.strip() if c.isalpha()]).upper()
        if name_clean:
            if name_clean[:2] in VDW_RADII:
                return VDW_RADII[name_clean[:2]]
            if name_clean[0] in VDW_RADII:
                return VDW_RADII[name_clean[0]]
        return DEFAULT_VDW_RADIUS


@dataclass
class Residue:
    chain: str
    resseq: int
    icode: str
    resname: str
    std_resname: str  # 標準20アミノ酸名（非タンパク質の場合は resname と同じ）
    is_protein: bool
    atoms: List[Atom] = field(default_factory=list)

    @property
    def key(self) -> Tuple[str, int, str]:
        return (self.chain, self.resseq, self.icode)

    @property
    def label(self) -> str:
        ic = self.icode.strip()
        return f"{self.resname}{self.resseq}{ic}"

    @property
    def heavy_atoms(self) -> List[Atom]:
        heavies = [a for a in self.atoms if not a.is_hydrogen]
        return heavies if heavies else self.atoms

    @property
    def ca_atom(self) -> Atom:
        for a in self.atoms:
            if a.name.strip() == "CA":
                return a
        return self.heavy_atoms[0]

    @property
    def centroid(self) -> Tuple[float, float, float]:
        pts = self.heavy_atoms
        n = len(pts)
        return (
            sum(a.x for a in pts) / n,
            sum(a.y for a in pts) / n,
            sum(a.z for a in pts) / n,
        )

    @property
    def sasa(self) -> float:
        return sum(a.sasa for a in self.atoms)

    @property
    def mean_bfactor(self) -> float:
        pts = self.heavy_atoms
        return sum(a.bfactor for a in pts) / len(pts) if pts else 0.0

    def rsa_tien2013(self, empirical: bool = False, cap_at_one: bool = True) -> float:
        scale = TIEN_2013_MAX_ASA_EMPIRICAL if empirical else TIEN_2013_MAX_ASA_THEORETICAL
        max_asa = scale.get(self.std_resname, 200.0)
        val = self.sasa / max_asa if max_asa > 0 else 0.0
        return min(1.0, val) if cap_at_one else val


@dataclass
class LigandEntity:
    role: str           # "Cofactor" or "Substrate" or "Ligand"
    resname: str
    chain: str
    resseq: int
    icode: str
    atoms: List[Atom]

    @property
    def label(self) -> str:
        ic = self.icode.strip()
        ch = self.chain.strip() or "_"
        return f"{self.resname}({ch}:{self.resseq}{ic})"

    @property
    def heavy_atoms(self) -> List[Atom]:
        heavies = [a for a in self.atoms if not a.is_hydrogen]
        return heavies if heavies else self.atoms

    @property
    def centroid(self) -> Tuple[float, float, float]:
        pts = self.heavy_atoms
        n = len(pts)
        return (
            sum(a.x for a in pts) / n,
            sum(a.y for a in pts) / n,
            sum(a.z for a in pts) / n,
        )


@dataclass
class TerminusMetrics:
    terminus_type: str        # "N-terminus" or "C-terminus"
    chain: str
    residue_label: str        # e.g., "MET1"
    resseq: int
    sasa_1res: float          # 末端残基の絶対 SASA (Å^2)
    rsa_1res: float           # 末端1残基の相対溶媒露出度 RSA (Tien et al., 2013)
    sasa_win3: float          # 末端3残基窓の平均 SASA (Å^2)
    rsa_win3: float           # 末端3残基窓の平均 RSA (Thornton & Sibanda, 1983; Tien et al., 2013)
    exposure_class: str       # "Exposed (露出)", "Partially Exposed (部分露出)", "Buried (埋没)"
    mean_bfactor: float       # 末端3残基の平均 B-factor / pLDDT
    dist_min_all_ligands: Optional[float]      # 全リガンドに対する最小重原子間距離 (Å)
    dist_ca_centroid_all: Optional[float]      # 全リガンド重心に対する末端 Cα 距離 (Å)
    nearest_ligand_label: str                  # 最も近いリガンドの識別子
    dist_min_substrate: Optional[float]        # 基質 (Substrate) との最小重原子間距離 (Å)
    dist_ca_substrate: Optional[float]         # 基質重心との Cα 距離 (Å)
    substrate_label: str
    dist_min_cofactor: Optional[float]         # 補酵素 (Cofactor) との最小重原子間距離 (Å)
    dist_ca_cofactor: Optional[float]          # 補酵素重心との Cα 距離 (Å)
    cofactor_label: str
    per_ligand_distances: Dict[str, Tuple[float, float]]  # label -> (d_min, d_ca_centroid)
    orientation_angle_deg: Optional[float]     # 末端進行ベクトルと活性部位方向ベクトルのなす角 (度)
    dist_risk_score: float    # 距離リスクスコア (0-100)
    rsa_risk_score: float     # 溶媒埋没リスクスコア (0-100)
    orient_risk_score: float  # 配向リスクスコア (0-100)
    total_risk_score: float   # 総合リスクスコア (0-100)
    risk_level: str           # "LOW (低リスク)", "MODERATE (中リスク)", "HIGH (高リスク)"
    design_comment: str       # 具体的なタグ・リンカー設計指針


@dataclass
class PDBEvaluationResult:
    pdb_name: str
    pdb_path: Path
    num_protein_chains: int
    protein_chains: List[str]
    total_residues: int
    ligands: List[LigandEntity]
    n_term: TerminusMetrics
    c_term: TerminusMetrics
    all_chain_termini: List[Tuple[TerminusMetrics, TerminusMetrics]]
    recommended_terminus: str
    recommendation_reason: str
    image_paths: List[Path]   # [view1_0deg.png, view2_rot90_1.png, view3_rot90_2.png]
    pml_path: Optional[Path] = None


# =============================================================================
# 3. PDB パーサー & 基質・補酵素自動分類 (PDB Parser & Ligand Classifier)
# =============================================================================

def parse_pdb_file(
    pdb_path: Path,
    substrate_resnames: Optional[Sequence[str]] = None,
    cofactor_resnames: Optional[Sequence[str]] = None,
    extra_ignore_resnames: Optional[Sequence[str]] = None,
) -> Tuple[Dict[str, List[Residue]], List[LigandEntity]]:
    """
    PDBファイルを解析し、タンパク質鎖ごとの残基リストと、リガンド（基質・補酵素）リストを返す。
    MODEL が複数ある場合は最初の MODEL のみを使用する。
    """
    sub_set = {r.strip().upper() for r in (substrate_resnames or []) if r.strip()}
    cof_set = {r.strip().upper() for r in (cofactor_resnames or []) if r.strip()}
    ignore_set = WATER_RESNAMES | BUFFER_AND_SALT_RESNAMES
    if extra_ignore_resnames:
        ignore_set |= {r.strip().upper() for r in extra_ignore_resnames if r.strip()}
    # 明示指定された基質・補酵素は無視リストから外す
    ignore_set -= (sub_set | cof_set)

    raw_atoms: List[Atom] = []
    in_first_model = True

    with open(pdb_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            rec = line[:6].strip().upper()
            if rec == "ENDMDL":
                if raw_atoms:
                    break
                in_first_model = False
                continue
            if not in_first_model:
                continue
            if rec not in ("ATOM", "HETATM"):
                continue
            if len(line) < 54:
                continue

            altloc = line[16:17].strip()
            # 代替コンフォメーションがある場合は A または 1 または空白のみ採用
            if altloc not in ("", "A", "1"):
                continue

            try:
                serial = int(line[6:11].strip() or "0")
            except ValueError:
                serial = 0
            name = line[12:16]
            resname = line[17:20].strip().upper()
            chain = line[21:22].strip() or "A"
            try:
                resseq = int(line[22:26].strip() or "0")
            except ValueError:
                resseq = 0
            icode = line[26:27].strip()
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            try:
                occupancy = float(line[54:60].strip()) if len(line) >= 60 and line[54:60].strip() else 1.0
            except ValueError:
                occupancy = 1.0
            try:
                bfactor = float(line[60:66].strip()) if len(line) >= 66 and line[60:66].strip() else 0.0
            except ValueError:
                bfactor = 0.0

            element = line[76:78].strip() if len(line) >= 78 else ""
            if not element:
                # 原子名から推定
                alpha = "".join([c for c in name.strip() if c.isalpha()])
                element = alpha[:1].upper() if alpha else "C"

            raw_atoms.append(
                Atom(
                    serial=serial,
                    name=name,
                    altloc=altloc,
                    resname=resname,
                    chain=chain,
                    resseq=resseq,
                    icode=icode,
                    x=x,
                    y=y,
                    z=z,
                    occupancy=occupancy,
                    bfactor=bfactor,
                    element=element,
                    record_type=rec,
                )
            )

    # 残基単位にグループ化
    residues_by_key: Dict[Tuple[str, int, str, str], List[Atom]] = {}
    residue_order: List[Tuple[str, int, str, str]] = []
    for atom in raw_atoms:
        # 同一残基番号でもタンパク質とHETATMリガンドが重ならないよう record_type + resname も考慮
        k = (atom.chain, atom.resseq, atom.icode, atom.resname)
        if k not in residues_by_key:
            residues_by_key[k] = []
            residue_order.append(k)
        residues_by_key[k].append(atom)

    protein_chains: Dict[str, List[Residue]] = {}
    candidate_ligands: List[Tuple[str, str, int, str, List[Atom]]] = []

    for chain, resseq, icode, resname in residue_order:
        atoms = residues_by_key[(chain, resseq, icode, resname)]
        # 明示的に基質・補酵素として指定されている場合はリガンド扱い
        if resname in sub_set or resname in cof_set:
            candidate_ligands.append((resname, chain, resseq, icode, atoms))
            continue

        # タンパク質残基判定 (ATOM レコードかつ標準アミノ酸、またはペプチド主鎖 N, CA, C を持つ標準アミノ酸)
        is_std_aa = resname in PROTEIN_RESIDUE_MAP
        atom_names = {a.name.strip().upper() for a in atoms}
        has_backbone = {"N", "CA", "C"}.issubset(atom_names)

        if is_std_aa and (atoms[0].record_type == "ATOM" or has_backbone):
            res = Residue(
                chain=chain,
                resseq=resseq,
                icode=icode,
                resname=resname,
                std_resname=PROTEIN_RESIDUE_MAP[resname],
                is_protein=True,
                atoms=atoms,
            )
            protein_chains.setdefault(chain, []).append(res)
        else:
            # 水や緩衝液イオンは除外
            if resname in ignore_set:
                continue
            # 重原子数が1個のみの単原子イオンで、かつ補酵素指定がない場合は除外（ただし触媒金属を明示した場合は残る）
            heavies = [a for a in atoms if not a.is_hydrogen]
            if len(heavies) <= 1 and resname not in cof_set and resname not in sub_set:
                continue
            candidate_ligands.append((resname, chain, resseq, icode, atoms))

    # 各タンパク質鎖の残基を resseq 順に整列（PDB記載順を基本としつつ確認）
    # リガンドの役割 (Substrate / Cofactor) の自動判別
    ligands: List[LigandEntity] = []
    for resname, chain, resseq, icode, atoms in candidate_ligands:
        if resname in sub_set:
            role = "Substrate"
        elif resname in cof_set or resname in KNOWN_COFACTOR_RESNAMES:
            role = "Cofactor"
        else:
            role = "Substrate"  # デフォルトはドッキング基質候補とする
        ligands.append(
            LigandEntity(
                role=role,
                resname=resname,
                chain=chain,
                resseq=resseq,
                icode=icode,
                atoms=atoms,
            )
        )

    # もし2種類以上のリガンドがあり、すべてが "Substrate" になった場合（両方とも UNL / LIG などの場合）、
    # 重原子数が多い方（または既知補酵素に近い方）と基質を区別できるようにしつつ両方評価対象にする
    if len(ligands) >= 2:
        roles = [lig.role for lig in ligands]
        if all(r == "Substrate" for r in roles):
            # 1つ目の非タンパク質リガンドと2つ目を区別しやすいよう、名称や重原子数で整理
            # (例: 補酵素が先に記載され、基質が後に記載されるケースやその逆に対応)
            # どちらも基質/補酵素候補として個別距離も全体最小距離も計算する
            ligands[0].role = "Substrate/Ligand1"
            ligands[1].role = "Cofactor/Ligand2"

    return protein_chains, ligands


# =============================================================================
# 4. Shrake-Rupley SASA & Tien et al. (2013) RSA 計算アルゴリズム
# =============================================================================

def _generate_fibonacci_sphere(n_points: int = 960) -> List[Tuple[float, float, float]]:
    """黄金螺旋 (Fibonacci sphere) により球面上に均一分布する単位ベクトル群を生成する。"""
    pts: List[Tuple[float, float, float]] = []
    phi = math.pi * (3.0 - math.sqrt(5.0))  # 黄金角
    for i in range(n_points):
        y = 1.0 - (i / float(n_points - 1)) * 2.0
        radius = math.sqrt(max(0.0, 1.0 - y * y))
        theta = phi * i
        x = math.cos(theta) * radius
        z = math.sin(theta) * radius
        pts.append((x, y, z))
    return pts


def compute_sasa_shrake_rupley(
    protein_chains: Dict[str, List[Residue]],
    n_sphere_points: int = 960,
    probe_radius: float = WATER_PROBE_RADIUS,
    target_residue_keys: Optional[set] = None,
) -> None:
    """
    Shrake & Rupley (1973) J. Mol. Biol. 79:351-371 アルゴリズムに基づき、
    アポ酵素（タンパク質重原子群）における各原子の溶媒露出表面積 (SASA, Å^2) を計算し、
    各 Atom.sasa に格納する。

    Parameters
    ----------
    protein_chains : Dict[str, List[Residue]]
        タンパク質鎖の辞書
    n_sphere_points : int
        各原子球面上のテスト点数 (デフォルト: 960)
    probe_radius : float
        水分子プローブ半径 (デフォルト: 1.40 Å)
    target_residue_keys : Optional[set]
        計算対象の残基キー集合。指定時はその残基に属する原子のみ球面テストを実行し、
        周囲の遮蔽原子には全タンパク質重原子を用いる（高速かつ厳密）。
    """
    all_heavy_atoms: List[Atom] = []
    target_atoms: List[Atom] = []

    for chain_residues in protein_chains.values():
        for res in chain_residues:
            for atom in res.atoms:
                atom.sasa = 0.0
                if atom.is_hydrogen:
                    continue
                all_heavy_atoms.append(atom)
                if target_residue_keys is None or res.key in target_residue_keys:
                    target_atoms.append(atom)

    if not all_heavy_atoms or not target_atoms:
        return

    sphere_pts = _generate_fibonacci_sphere(n_sphere_points)
    four_pi_over_n = (4.0 * math.pi) / float(n_sphere_points)

    # 空間ハッシュ（ボクセルグリッド）を構築して近傍探索を O(1) 化
    max_expanded_r = max(a.vdw_radius for a in all_heavy_atoms) + probe_radius
    cell_size = max_expanded_r * 2.0

    grid: Dict[Tuple[int, int, int], List[Tuple[Atom, float, float]]] = {}
    for atom in all_heavy_atoms:
        r_exp = atom.vdw_radius + probe_radius
        r_exp_sq = r_exp * r_exp
        gx = int(math.floor(atom.x / cell_size))
        gy = int(math.floor(atom.y / cell_size))
        gz = int(math.floor(atom.z / cell_size))
        grid.setdefault((gx, gy, gz), []).append((atom, r_exp, r_exp_sq))

    for atom_i in target_atoms:
        ri = atom_i.vdw_radius + probe_radius
        xi, yi, zi = atom_i.x, atom_i.y, atom_i.z
        gx = int(math.floor(xi / cell_size))
        gy = int(math.floor(yi / cell_size))
        gz = int(math.floor(zi / cell_size))

        # 近傍原子の収集
        neighbors: List[Tuple[float, float, float, float]] = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    cell = grid.get((gx + dx, gy + dy, gz + dz))
                    if not cell:
                        continue
                    for atom_j, rj, rj_sq in cell:
                        if atom_j is atom_i:
                            continue
                        ddx = atom_j.x - xi
                        ddy = atom_j.y - yi
                        ddz = atom_j.z - zi
                        dist_sq = ddx * ddx + ddy * ddy + ddz * ddz
                        max_d = ri + rj
                        if dist_sq < max_d * max_d:
                            neighbors.append((atom_j.x, atom_j.y, atom_j.z, rj_sq))

        if not neighbors:
            atom_i.sasa = 4.0 * math.pi * ri * ri
            continue

        # 近い順にソートすると早期棄却が高速化される
        neighbors.sort(key=lambda nb: (nb[0] - xi) ** 2 + (nb[1] - yi) ** 2 + (nb[2] - zi) ** 2)

        accessible_count = 0
        for sx, sy, sz in sphere_pts:
            px = xi + ri * sx
            py = yi + ri * sy
            pz = zi + ri * sz
            buried = False
            for nx, ny, nz, rj_sq in neighbors:
                dx = px - nx
                dy = py - ny
                dz = pz - nz
                if dx * dx + dy * dy + dz * dz < rj_sq:
                    buried = True
                    break
            if not buried:
                accessible_count += 1

        atom_i.sasa = four_pi_over_n * (ri * ri) * accessible_count


# =============================================================================
# 5. 距離・配向角・文献基準に基づくリスクスコア算出 (Risk Evaluation Engine)
# =============================================================================

def _min_heavy_atom_distance(res: Residue, lig_atoms: Sequence[Atom]) -> float:
    """残基の重原子群とリガンド重原子群の間の最小ユークリッド距離 (Å) を返す。"""
    r_atoms = res.heavy_atoms
    l_atoms = [a for a in lig_atoms if not a.is_hydrogen]
    if not l_atoms:
        l_atoms = list(lig_atoms)
    min_d_sq = float("inf")
    for ra in r_atoms:
        rx, ry, rz = ra.x, ra.y, ra.z
        for la in l_atoms:
            dx = rx - la.x
            dy = ry - la.y
            dz = rz - la.z
            d_sq = dx * dx + dy * dy + dz * dz
            if d_sq < min_d_sq:
                min_d_sq = d_sq
    return math.sqrt(min_d_sq)


def _point_distance(p1: Tuple[float, float, float], p2: Tuple[float, float, float]) -> float:
    return math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2 + (p1[2] - p2[2]) ** 2)


def _compute_terminal_orientation_angle(
    chain_residues: List[Residue],
    is_n_term: bool,
    active_site_centroid: Optional[Tuple[float, float, float]],
) -> Optional[float]:
    """
    末端から外側へ向かう主鎖ベクトル v_term と、末端から活性部位（リガンド重心）へ向かう
    ベクトル v_lig のなす角度 θ (0°〜180°) を計算する。
    - θ が小さい (< 60°): 末端の伸長方向が活性部位ポケットへ向いている（高リスク傾向）
    - θ が大きい (> 90°): 末端の伸長方向が活性部位と反対側の溶媒方向へ向いている（低リスク傾向）
    """
    if active_site_centroid is None or len(chain_residues) < 2:
        return None

    if is_n_term:
        r_term = chain_residues[0].ca_atom
        r_inner = chain_residues[min(2, len(chain_residues) - 1)].ca_atom
    else:
        r_term = chain_residues[-1].ca_atom
        r_inner = chain_residues[max(0, len(chain_residues) - 3)].ca_atom

    # 内側から末端へ向かうベクトル (タグが付加された際に伸びていく初期方向)
    vx = r_term.x - r_inner.x
    vy = r_term.y - r_inner.y
    vz = r_term.z - r_inner.z
    v_norm = math.sqrt(vx * vx + vy * vy + vz * vz)
    if v_norm < 1e-6:
        return None

    # 末端からリガンド重心へ向かうベクトル
    lx = active_site_centroid[0] - r_term.x
    ly = active_site_centroid[1] - r_term.y
    lz = active_site_centroid[2] - r_term.z
    l_norm = math.sqrt(lx * lx + ly * ly + lz * lz)
    if l_norm < 1e-6:
        return 0.0

    cos_theta = max(-1.0, min(1.0, (vx * lx + vy * ly + vz * lz) / (v_norm * l_norm)))
    return math.degrees(math.acos(cos_theta))


def score_distance_risk(d_min: Optional[float]) -> float:
    """
    リガンド（基質・補酵素）との最小重原子間距離 d_min (Å) から距離リスクスコア (0–100) を算出。
    根拠文献:
      - Carson et al. (2007) Acta Crystallogr. D 63:295-301:
        His6タグの完全伸長時の長さは約 3.6 Å × 6 ≈ 21.6 Å (ランダムコイル到達半径 ~12–15 Å)。
      - Freydank et al. (2008) Proteins 72:173-183 & Yeon et al. (2014) Biotechnol. Bioprocess Eng. 19:798-802:
        末端が活性部位から 15 Å 未満にある場合、His6タグによる基質結合阻害・触媒ループ変形の危険が高い。
        15–22 Å は最大伸長圏内（中リスク）、22 Å 以上は到達圏外（低リスク）。
    """
    if d_min is None:
        return 25.0  # リガンド座標なしの場合は中立値
    if d_min <= 8.0:
        return 100.0
    if d_min < 15.0:
        # 8 Å -> 100点, 15 Å -> 70点 (高リスク領域)
        return 70.0 + (15.0 - d_min) / (15.0 - 8.0) * 30.0
    if d_min < 22.0:
        # 15 Å -> 70点, 22 Å -> 30点 (中リスク領域: His6最大伸長長圏内)
        return 30.0 + (22.0 - d_min) / (22.0 - 15.0) * 40.0
    if d_min < 32.0:
        # 22 Å -> 30点, 32 Å -> 0点 (低リスク領域)
        return max(0.0, (32.0 - d_min) / (32.0 - 22.0) * 30.0)
    return 0.0


def score_rsa_risk(rsa_win3: float) -> Tuple[float, str]:
    """
    末端領域の相対溶媒露出度 RSA から溶媒埋没リスクスコア (0–100) と露出度クラスを算出。
    根拠文献:
      - Tien et al. (2013) PLOS ONE 8(11): e80635 (MaxASA基準値)
      - Miller et al. (1987) J. Mol. Biol. 196:641-656 & Rost & Sander (1994) Proteins 20:216-226
      - Thornton & Sibanda (1983) J. Mol. Biol. 167:443-460 (末端領域の露出性とフォールディング安定性)
        * RSA < 0.20 (20%未満): Buried (埋没) -> タグ付加で立体構造コアの崩壊・封入体形成・Ni-NTA非結合リスク大
        * 0.20 <= RSA < 0.50: Partially Exposed (部分露出) -> フレキシブルリンカー (GGGS等) 挿入が望ましい
        * RSA >= 0.50 (50%以上): Exposed (溶媒露出) -> タグ付加に好適
    """
    if rsa_win3 < 0.20:
        exposure_class = "Buried (埋没: <20%)"
        # RSA 0.0 -> 100点, RSA 0.20 -> 70点
        score = 70.0 + (0.20 - rsa_win3) / 0.20 * 30.0
    elif rsa_win3 < 0.50:
        exposure_class = "Partially Exposed (部分露出: 20-50%)"
        # RSA 0.20 -> 70点, RSA 0.50 -> 25点
        score = 25.0 + (0.50 - rsa_win3) / (0.50 - 0.20) * 45.0
    else:
        exposure_class = "Exposed (露出: >=50%)"
        # RSA 0.50 -> 25点, RSA >= 0.85 -> 0点
        score = max(0.0, (0.85 - min(0.85, rsa_win3)) / (0.85 - 0.50) * 25.0)
    return round(score, 1), exposure_class


def score_orientation_risk(angle_deg: Optional[float]) -> float:
    """
    末端ベクトルと活性部位方向ベクトルのなす角 θ (0–180°) から配向リスクスコア (0–100) を算出。
    - θ <= 45°: 活性部位の方向へ直接向かっている -> 80–100点
    - θ >= 120°: 活性部位と反対方向へ向かっている -> 0–15点
    """
    if angle_deg is None:
        return 25.0
    return round(max(0.0, min(100.0, (1.0 - angle_deg / 150.0) * 100.0)), 1)


def evaluate_terminus(
    chain_residues: List[Residue],
    is_n_term: bool,
    ligands: List[LigandEntity],
) -> TerminusMetrics:
    """1本の鎖の N末端 または C末端 の各種構造指標と Hisタグ付加リスクを評価する。"""
    term_type = "N-terminus" if is_n_term else "C-terminus"
    term_res = chain_residues[0] if is_n_term else chain_residues[-1]
    win_residues = chain_residues[:3] if is_n_term else chain_residues[-3:]

    sasa_1 = term_res.sasa
    rsa_1 = term_res.rsa_tien2013()
    sasa_w3 = sum(r.sasa for r in win_residues) / len(win_residues)
    rsa_w3 = sum(r.rsa_tien2013() for r in win_residues) / len(win_residues)
    mean_b = sum(r.mean_bfactor for r in win_residues) / len(win_residues)

    # リガンドとの距離計算
    per_lig: Dict[str, Tuple[float, float]] = {}
    d_min_all: Optional[float] = None
    d_ca_all: Optional[float] = None
    nearest_label = "None"

    d_min_sub: Optional[float] = None
    d_ca_sub: Optional[float] = None
    sub_label = "None"

    d_min_cof: Optional[float] = None
    d_ca_cof: Optional[float] = None
    cof_label = "None"

    term_ca_pt = (term_res.ca_atom.x, term_res.ca_atom.y, term_res.ca_atom.z)

    all_lig_atoms: List[Atom] = []
    for lig in ligands:
        all_lig_atoms.extend(lig.heavy_atoms)
        d_m = _min_heavy_atom_distance(term_res, lig.heavy_atoms)
        d_c = _point_distance(term_ca_pt, lig.centroid)
        per_lig[lig.label] = (round(d_m, 2), round(d_c, 2))

        if d_min_all is None or d_m < d_min_all:
            d_min_all = d_m
            d_ca_all = d_c
            nearest_label = f"{lig.role}:{lig.label}"

        if "Cofactor" in lig.role:
            if d_min_cof is None or d_m < d_min_cof:
                d_min_cof = d_m
                d_ca_cof = d_c
                cof_label = lig.label
        else:
            if d_min_sub is None or d_m < d_min_sub:
                d_min_sub = d_m
                d_ca_sub = d_c
                sub_label = lig.label

    # 活性部位全体の重心
    if all_lig_atoms:
        active_centroid: Optional[Tuple[float, float, float]] = (
            sum(a.x for a in all_lig_atoms) / len(all_lig_atoms),
            sum(a.y for a in all_lig_atoms) / len(all_lig_atoms),
            sum(a.z for a in all_lig_atoms) / len(all_lig_atoms),
        )
    else:
        active_centroid = None

    angle_deg = _compute_terminal_orientation_angle(chain_residues, is_n_term, active_centroid)

    # スコア算出
    # 重み付け: 活性部位距離リスク 55% + 溶媒埋没リスク 35% + 配向リスク 10%
    # ただし距離が極端に近い (< 12 Å) または完全に埋没 (RSA < 0.12) の場合はペナルティで HIGH リスクを保証
    d_score = score_distance_risk(d_min_all)
    r_score, exp_class = score_rsa_risk(rsa_w3)
    o_score = score_orientation_risk(angle_deg)

    weighted = 0.55 * d_score + 0.35 * r_score + 0.10 * o_score
    # ボトルネック補正（距離または埋没のどちらかが致命的なら総合リスクも引き上げる）
    bottleneck = max(d_score, r_score) * 0.85
    total_score = round(min(100.0, max(weighted, bottleneck)), 1)

    if total_score >= 60.0 or (d_min_all is not None and d_min_all < 15.0) or rsa_w3 < 0.20:
        risk_level = "HIGH (高リスク)"
    elif total_score >= 32.0 or (d_min_all is not None and d_min_all < 22.0) or rsa_w3 < 0.40:
        risk_level = "MODERATE (中リスク)"
    else:
        risk_level = "LOW (低リスク)"

    # 設計コメントの生成
    reasons: List[str] = []
    if d_min_all is not None:
        if d_min_all < 15.0:
            reasons.append(f"リガンドまで {d_min_all:.1f} Å (<15 Å: His6直接干渉圏内)")
        elif d_min_all < 22.0:
            reasons.append(f"リガンドまで {d_min_all:.1f} Å (15-22 Å: His6最大伸長圏内)")
        else:
            reasons.append(f"リガンドまで十分離れている ({d_min_all:.1f} Å >= 22 Å)")
    if rsa_w3 < 0.20:
        reasons.append(f"末端が埋没 (RSA={rsa_w3*100:.1f}% < 20%: 構造不安定化/精製不良リスク)")
    elif rsa_w3 < 0.50:
        reasons.append(f"末端が部分露出 (RSA={rsa_w3*100:.1f}%: GGGS等の短いリンカー推奨)")
    else:
        reasons.append(f"末端が十分溶媒露出 (RSA={rsa_w3*100:.1f}% >= 50%)")

    return TerminusMetrics(
        terminus_type=term_type,
        chain=term_res.chain,
        residue_label=term_res.label,
        resseq=term_res.resseq,
        sasa_1res=round(sasa_1, 2),
        rsa_1res=round(rsa_1, 3),
        sasa_win3=round(sasa_w3, 2),
        rsa_win3=round(rsa_w3, 3),
        exposure_class=exp_class,
        mean_bfactor=round(mean_b, 2),
        dist_min_all_ligands=round(d_min_all, 2) if d_min_all is not None else None,
        dist_ca_centroid_all=round(d_ca_all, 2) if d_ca_all is not None else None,
        nearest_ligand_label=nearest_label,
        dist_min_substrate=round(d_min_sub, 2) if d_min_sub is not None else None,
        dist_ca_substrate=round(d_ca_sub, 2) if d_ca_sub is not None else None,
        substrate_label=sub_label,
        dist_min_cofactor=round(d_min_cof, 2) if d_min_cof is not None else None,
        dist_ca_cofactor=round(d_ca_cof, 2) if d_ca_cof is not None else None,
        cofactor_label=cof_label,
        per_ligand_distances=per_lig,
        orientation_angle_deg=round(angle_deg, 1) if angle_deg is not None else None,
        dist_risk_score=round(d_score, 1),
        rsa_risk_score=round(r_score, 1),
        orient_risk_score=round(o_score, 1),
        total_risk_score=total_score,
        risk_level=risk_level,
        design_comment=" / ".join(reasons),
    )


def evaluate_pdb_structure(
    pdb_path: Path,
    substrate_resnames: Optional[Sequence[str]] = None,
    cofactor_resnames: Optional[Sequence[str]] = None,
) -> PDBEvaluationResult:
    """単一のPDBファイルを解析し、SASA/RSA・距離・Hisタグ付加リスク評価結果を返す。"""
    protein_chains, ligands = parse_pdb_file(
        pdb_path,
        substrate_resnames=substrate_resnames,
        cofactor_resnames=cofactor_resnames,
    )
    if not protein_chains:
        raise ValueError(f"PDBファイル内にタンパク質鎖 (ATOMレコード) が見つかりません: {pdb_path}")

    # 末端3残基のキーを収集して Shrake-Rupley SASA を高速・厳密に計算
    target_keys = set()
    for residues in protein_chains.values():
        for r in residues[:3] + residues[-3:]:
            target_keys.add(r.key)

    compute_sasa_shrake_rupley(protein_chains, n_sphere_points=960, target_residue_keys=target_keys)

    all_chain_termini: List[Tuple[TerminusMetrics, TerminusMetrics]] = []
    for chain_id, residues in protein_chains.items():
        n_m = evaluate_terminus(residues, is_n_term=True, ligands=ligands)
        c_m = evaluate_terminus(residues, is_n_term=False, ligands=ligands)
        all_chain_termini.append((n_m, c_m))

    # 代表鎖（最初の酵素鎖、通常 Chain A）の N末・C末
    primary_n, primary_c = all_chain_termini[0]

    # 推奨タグ位置の決定
    diff = primary_c.total_risk_score - primary_n.total_risk_score
    if primary_n.risk_level.startswith("HIGH") and primary_c.risk_level.startswith("HIGH"):
        better = "N末" if primary_n.total_risk_score <= primary_c.total_risk_score else "C末"
        rec = f"両末端とも高リスク ({better}側 + 切断可能リンカー推奨)"
        reason = (
            f"N末(Score={primary_n.total_risk_score})・C末(Score={primary_c.total_risk_score})ともにリスクが高いため、"
            f"比較的スコアの低い{better}側にフレキシブルリンカー(GGGS)2およびTEV等のプロテアーゼ切断配列を介してHisタグを付加することを推奨。"
        )
    elif primary_n.total_risk_score + 5.0 < primary_c.total_risk_score:
        rec = "N-terminus (N末推奨)"
        reason = (
            f"N末(Score={primary_n.total_risk_score}, {primary_n.risk_level}) の方が "
            f"C末(Score={primary_c.total_risk_score}, {primary_c.risk_level}) よりも低リスク。"
        )
    elif primary_c.total_risk_score + 5.0 < primary_n.total_risk_score:
        rec = "C-terminus (C末推奨)"
        reason = (
            f"C末(Score={primary_c.total_risk_score}, {primary_c.risk_level}) の方が "
            f"N末(Score={primary_n.total_risk_score}, {primary_n.risk_level}) よりも低リスク。"
        )
    else:
        better = "C-terminus (C末)" if primary_c.total_risk_score <= primary_n.total_risk_score else "N-terminus (N末)"
        rec = f"Both Comparable ({better} やや優位)"
        reason = (
            f"N末(Score={primary_n.total_risk_score}) と C末(Score={primary_c.total_risk_score}) は同程度のリスク水準。"
        )

    total_res = sum(len(v) for v in protein_chains.values())
    return PDBEvaluationResult(
        pdb_name=pdb_path.stem,
        pdb_path=pdb_path,
        num_protein_chains=len(protein_chains),
        protein_chains=list(protein_chains.keys()),
        total_residues=total_res,
        ligands=ligands,
        n_term=primary_n,
        c_term=primary_c,
        all_chain_termini=all_chain_termini,
        recommended_terminus=rec,
        recommendation_reason=reason,
        image_paths=[],
    )


# =============================================================================
# 6. Windows デスクトップ版 PyMOL / WSL2 Ubuntu ブリッジ & 3方向バッチ作図
# =============================================================================

def is_wsl2() -> bool:
    """現在のPython実行環境が WSL / WSL2 (Ubuntu等) かどうかを判定する。"""
    if os.environ.get("WSL_DISTRO_NAME") or os.environ.get("WSL_INTEROP"):
        return True
    try:
        with open("/proc/version", "r", encoding="utf-8") as f:
            ver = f.read().lower()
            return "microsoft" in ver or "wsl" in ver
    except OSError:
        return False


def find_pymol_executable(user_specified_cmd: Optional[str] = None) -> Optional[str]:
    """
    Windows 環境または WSL2 (Ubuntu) 環境から、利用可能な PyMOL 実行ファイルを自動検出する。
    Windows PC にインストールされたデスクトップ版 PyMOL (Schrodinger PyMOL2 / PyMOLWin.exe / pymol.exe)
    を WSL2 からもシームレスに呼び出せるように探索する。
    """
    if user_specified_cmd:
        if os.path.isfile(user_specified_cmd) or shutil.which(user_specified_cmd):
            return user_specified_cmd

    env_pymol = os.environ.get("PYMOL_PATH") or os.environ.get("PYMOL_EXE")
    if env_pymol and (os.path.isfile(env_pymol) or shutil.which(env_pymol)):
        return env_pymol

    # PATH 上の pymol / PyMOLWin.exe
    for cmd_name in ("pymol", "PyMOLWin.exe", "pymol.exe"):
        found = shutil.which(cmd_name)
        if found:
            return found

    candidates: List[str] = []
    if platform.system() == "Windows":
        local_app = os.environ.get("LOCALAPPDATA", "")
        prog_data = os.environ.get("PROGRAMDATA", r"C:\ProgramData")
        prog_files = os.environ.get("ProgramFiles", r"C:\Program Files")
        prog_files_x86 = os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)")
        user_prof = os.environ.get("USERPROFILE", "")

        search_patterns = [
            os.path.join(local_app, "Schrodinger", "PyMOL2", "PyMOLWin.exe"),
            os.path.join(local_app, "Schrodinger", "PyMOL2", "pymol.exe"),
            os.path.join(prog_data, "Schrodinger", "PyMOL2", "PyMOLWin.exe"),
            os.path.join(prog_data, "Schrodinger", "PyMOL2", "pymol.exe"),
            os.path.join(prog_files, "Schrodinger", "PyMOL2", "PyMOLWin.exe"),
            os.path.join(prog_files, "PyMOL", "PyMOLWin.exe"),
            os.path.join(prog_files_x86, "PyMOL", "PyMOLWin.exe"),
            os.path.join(user_prof, "Schrodinger", "PyMOL2", "PyMOLWin.exe"),
            os.path.join(user_prof, "anaconda3", "Library", "bin", "pymol.exe"),
            os.path.join(user_prof, "miniconda3", "Library", "bin", "pymol.exe"),
            os.path.join(user_prof, "AppData", "Local", "anaconda3", "Library", "bin", "pymol.exe"),
            os.path.join(user_prof, "AppData", "Local", "miniconda3", "Library", "bin", "pymol.exe"),
        ]
        candidates.extend(search_patterns)
    else:
        # WSL2 Ubuntu から Windows 側のデスクトップ版 PyMOL を探索 (/mnt/c/...)
        wsl_patterns = [
            "/mnt/c/Users/*/AppData/Local/Schrodinger/PyMOL2/PyMOLWin.exe",
            "/mnt/c/Users/*/AppData/Local/Schrodinger/PyMOL2/pymol.exe",
            "/mnt/c/ProgramData/Schrodinger/PyMOL2/PyMOLWin.exe",
            "/mnt/c/ProgramData/Schrodinger/PyMOL2/pymol.exe",
            "/mnt/c/Program Files/Schrodinger/PyMOL2/PyMOLWin.exe",
            "/mnt/c/Program Files/PyMOL/PyMOLWin.exe",
            "/mnt/c/Program Files (x86)/PyMOL/PyMOLWin.exe",
            "/mnt/c/Users/*/anaconda3/Library/bin/pymol.exe",
            "/mnt/c/Users/*/miniconda3/Library/bin/pymol.exe",
        ]
        for pat in wsl_patterns:
            candidates.extend(sorted(glob.glob(pat)))

    for cand in candidates:
        if cand and os.path.isfile(cand):
            return cand

    return None


def to_pymol_compatible_path(path: Path, pymol_exe: Optional[str]) -> str:
    """
    WSL2 から Windows 版 PyMOL (.exe) を呼び出す場合、Windows 側の PyMOL が解釈できる
    Windows パス (例: C:/Users/...) に自動変換する。
    """
    abs_path = path.resolve()
    is_win_exe = bool(pymol_exe and pymol_exe.lower().endswith(".exe"))

    if is_wsl2() and is_win_exe:
        try:
            res = subprocess.run(
                ["wslpath", "-m", str(abs_path)],
                capture_output=True,
                text=True,
                check=True,
            )
            win_path = res.stdout.strip()
            if win_path:
                return win_path
        except Exception:
            pass
        # /mnt/c/... の手動フォールバック変換
        s = str(abs_path)
        if s.startswith("/mnt/") and len(s) > 6 and s[6] == "/":
            drive = s[5].upper()
            return f"{drive}:{s[6:]}"

    return str(abs_path).replace("\\", "/")


def build_pymol_script_content(
    pdb_file_for_pymol: str,
    out_img1: str,
    out_img2: str,
    out_img3: str,
    result: PDBEvaluationResult,
    rotation_mode: str = "y90_x90",
    image_width: int = 1200,
    image_height: int = 900,
    ray_trace: bool = True,
) -> str:
    """
    指定された配色と回転条件に従う PyMOL スクリプト (.pml) を生成する。
      - 酵素 (タンパク質鎖): 青色 (blue)
      - リガンド (基質・補酵素の2種など): 黄色 (yellow)
      - N末端: 赤色 (red)
      - C末端: 緑色 (green)
      - 合計3種の図:
          View 1: 基準向き (orient)
          View 2: 90度回転 #1 (Y軸周り 90° 回転)
          View 3: 90度回転 #2 (X軸周り 90° 回転 または Y軸周りさらに 90° 回転)
    """
    # 各鎖の N末・C末セレクション式の作成
    n_sel_parts: List[str] = []
    c_sel_parts: List[str] = []
    n_cartoon_parts: List[str] = []
    c_cartoon_parts: List[str] = []

    for n_m, c_m in result.all_chain_termini:
        ch = n_m.chain
        n_sel_parts.append(f"(chain {ch} and resi {n_m.resseq})")
        c_sel_parts.append(f"(chain {ch} and resi {c_m.resseq})")
        # カートゥーン上でも末端が視認しやすいよう末端3残基を色付けし、最末端残基を球表示
        n_cartoon_parts.append(f"(chain {ch} and resi {n_m.resseq}-{n_m.resseq + 2})")
        c_cartoon_parts.append(f"(chain {ch} and resi {max(1, c_m.resseq - 2)}-{c_m.resseq})")

    n_term_expr = " or ".join(n_sel_parts) if n_sel_parts else "none"
    c_term_expr = " or ".join(c_sel_parts) if c_sel_parts else "none"
    n_region_expr = " or ".join(n_cartoon_parts) if n_cartoon_parts else n_term_expr
    c_region_expr = " or ".join(c_cartoon_parts) if c_cartoon_parts else c_term_expr

    # リガンド（基質・補酵素）のセレクション式
    lig_sel_parts: List[str] = []
    for lig in result.ligands:
        ch_part = f"chain {lig.chain} and " if lig.chain.strip() else ""
        lig_sel_parts.append(f"({ch_part}resn {lig.resname} and resi {lig.resseq})")
    if lig_sel_parts:
        lig_expr = " or ".join(lig_sel_parts)
    else:
        lig_expr = "organic or (hetatm and not solvent and not (resn HOH+WAT+NA+CL+K+SO4+PO4+GOL+EDO))"

    ray_flag = "1" if ray_trace else "0"

    if rotation_mode == "y90_y180":
        rot2_cmd = "turn y, 90"
        rot3_cmd = "turn y, 90"
    else:
        # デフォルト y90_x90: 1枚目を基準に Y軸90度回転 (側面図) と X軸90度回転 (上面/底面図) の直交3面図
        rot2_cmd = "turn y, 90"
        rot3_cmd = "turn y, -90\nturn x, 90"

    pml = f"""# Auto-generated PyMOL script for {result.pdb_name}
# Coloring: Enzyme = blue, Ligands (Substrate & Cofactor) = yellow, N-terminus = red, C-terminus = green
reinitialize
bg_color white
set ray_opaque_background, on
set antialias, 2
set orthoscopic, on
set depth_cue, 0
set specular, 0.25
set cartoon_fancy_helices, 1
set sphere_scale, 0.75
set stick_radius, 0.25

load "{pdb_file_for_pymol}", complex
hide everything, all

# 1. セレクション定義
select enzyme, polymer.protein
select ligands, ({lig_expr}) and not polymer.protein
select n_term_res, enzyme and ({n_term_expr})
select c_term_res, enzyme and ({c_term_expr})
select n_term_region, enzyme and ({n_region_expr})
select c_term_region, enzyme and ({c_region_expr})

# 2. 表示スタイル設定
show cartoon, enzyme
show sticks, ligands
show spheres, n_term_res and (name CA+N)
show sticks, n_term_res
show spheres, c_term_res and (name CA+C+OXT)
show sticks, c_term_res

# 3. カラーリング (酵素=青, リガンド[基質・補酵素]=黄, N末=赤, C末=緑)
color blue, enzyme
color yellow, ligands
color red, n_term_region
color red, n_term_res
color green, c_term_region
color green, c_term_res

deselect

# 4. 3方向の図 (0°正面, 90°回転1, 90°回転2) をバッチ出力
orient enzyme
zoom all, 3

# View 1: Initial Oriented View (0 deg)
png "{out_img1}", width={image_width}, height={image_height}, dpi=300, ray={ray_flag}

# View 2: Rotated 90 deg (#1)
{rot2_cmd}
png "{out_img2}", width={image_width}, height={image_height}, dpi=300, ray={ray_flag}

# View 3: Rotated 90 deg (#2)
{rot3_cmd}
png "{out_img3}", width={image_width}, height={image_height}, dpi=300, ray={ray_flag}

quit
"""
    return pml


def _render_fallback_3d_images(
    pdb_path: Path,
    result: PDBEvaluationResult,
    out_paths: List[Path],
    rotation_mode: str = "y90_x90",
) -> None:
    """
    PyMOL が未インストールの環境やヘッドレス検証環境（--fallback-render 等）向けに、
    同一の配色ルール（酵素=青、リガンド[基質・補酵素]=黄、N末=赤、C末=緑）かつ
    同一の90度回転（0°、90°回転#1、90°回転#2）を適用した3方向の構造投影図 (PNG) を生成する。
    """
    protein_chains, ligands = parse_pdb_file(pdb_path)

    # 3つの回転行列 (View 1: 0°, View 2: Y軸+90°, View 3: X軸+90° または Y軸+180°)
    def rot_coords(x: float, y: float, z: float, view_idx: int) -> Tuple[float, float, float]:
        if view_idx == 0:
            return (x, y, z)
        if view_idx == 1:
            # Y軸周り 90° 回転: (x, y, z) -> (z, y, -x)
            return (z, y, -x)
        if rotation_mode == "y90_y180":
            # Y軸周り 180° 回転: (x, y, z) -> (-x, y, -z)
            return (-x, y, -z)
        # X軸周り 90° 回転: (x, y, z) -> (x, -z, y)
        return (x, -z, y)

    # 全体の中心座標を計算して原点中心化
    all_ca = [res.ca_atom for residues in protein_chains.values() for res in residues]
    if not all_ca:
        return
    cx = sum(a.x for a in all_ca) / len(all_ca)
    cy = sum(a.y for a in all_ca) / len(all_ca)
    cz = sum(a.z for a in all_ca) / len(all_ca)

    view_titles = [
        f"{result.pdb_name} - View 1 (0 deg Front)",
        f"{result.pdb_name} - View 2 (Rot 90 deg Y)",
        f"{result.pdb_name} - View 3 (Rot 90 deg {'Y180' if rotation_mode == 'y90_y180' else 'X'})",
    ]

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        for v_idx, (out_p, title) in enumerate(zip(out_paths, view_titles)):
            fig = plt.figure(figsize=(5.2, 4.0), dpi=200, facecolor="white")
            ax = fig.add_subplot(111, projection="3d", facecolor="white")

            # 1. 酵素主鎖 (青色: blue)
            for ch_id, residues in protein_chains.items():
                xs, ys, zs = [], [], []
                for res in residues:
                    ca = res.ca_atom
                    rx, ry, rz = rot_coords(ca.x - cx, ca.y - cy, ca.z - cz, v_idx)
                    xs.append(rx)
                    ys.append(ry)
                    zs.append(rz)
                ax.plot(xs, ys, zs, color="#1565C0", linewidth=2.4, alpha=0.85, label="Enzyme (Blue)" if ch_id == result.protein_chains[0] else "")

            # 2. リガンド: 基質・補酵素 (黄色: yellow)
            for l_i, lig in enumerate(ligands):
                lxs, lys, lzs = [], [], []
                for a in lig.heavy_atoms:
                    rx, ry, rz = rot_coords(a.x - cx, a.y - cy, a.z - cz, v_idx)
                    lxs.append(rx)
                    lys.append(ry)
                    lzs.append(rz)
                ax.scatter(
                    lxs, lys, lzs,
                    color="#FFD600", edgecolors="#B28900", s=55, depthshade=True,
                    label=f"Ligand: {lig.role} ({lig.resname})" if l_i < 2 else "",
                    zorder=5,
                )

            # 3. N末端 (赤色: red) & C末端 (緑色: green)
            for c_i, (ch_id, residues) in enumerate(protein_chains.items()):
                n_ca = residues[0].ca_atom
                c_ca = residues[-1].ca_atom
                nx, ny, nz = rot_coords(n_ca.x - cx, n_ca.y - cy, n_ca.z - cz, v_idx)
                gx, gy, gz = rot_coords(c_ca.x - cx, c_ca.y - cy, c_ca.z - cz, v_idx)
                ax.scatter([nx], [ny], [nz], color="#D50000", edgecolors="black", s=130, zorder=10,
                           label=f"N-term ({residues[0].label})" if c_i == 0 else "")
                ax.scatter([gx], [gy], [gz], color="#00C853", edgecolors="black", s=130, zorder=10,
                           label=f"C-term ({residues[-1].label})" if c_i == 0 else "")

            ax.set_title(title, fontsize=9.5, fontweight="bold", pad=6)
            ax.legend(loc="upper right", fontsize=6.5, framealpha=0.85)
            ax.set_xlabel("X (Å)", fontsize=7)
            ax.set_ylabel("Y (Å)", fontsize=7)
            ax.set_zlabel("Z (Å)", fontsize=7)
            ax.tick_params(labelsize=6)
            plt.tight_layout()
            fig.savefig(str(out_p), dpi=200, facecolor="white")
            plt.close(fig)
        return
    except Exception:
        pass

    # Pure-Python PNG 2D/3D レンダラー (matplotlib すらない最小環境用の完全内蔵フォールバック)
    for v_idx, out_p in enumerate(out_paths):
        _write_pure_python_structure_png(
            out_p, protein_chains, ligands, (cx, cy, cz), v_idx, rotation_mode
        )


def _write_pure_python_structure_png(
    out_path: Path,
    protein_chains: Dict[str, List[Residue]],
    ligands: List[LigandEntity],
    center: Tuple[float, float, float],
    view_idx: int,
    rotation_mode: str,
    width: int = 480,
    height: int = 360,
) -> None:
    """外部ライブラリ非依存の Pure Python 3D投影 PNG 描画関数。"""
    cx, cy, cz = center

    def rot(x: float, y: float, z: float) -> Tuple[float, float, float]:
        if view_idx == 0:
            return (x, y, z)
        if view_idx == 1:
            return (z, y, -x)
        if rotation_mode == "y90_y180":
            return (-x, y, -z)
        return (x, -z, y)

    # 白背景 (RGB)
    buf = bytearray([255] * (width * height * 3))

    def set_px(px: int, py: int, r: int, g: int, b: int) -> None:
        if 0 <= px < width and 0 <= py < height:
            idx = (py * width + px) * 3
            buf[idx] = r
            buf[idx + 1] = g
            buf[idx + 2] = b

    def draw_circle(px: int, py: int, radius: int, rgb: Tuple[int, int, int], outline: Optional[Tuple[int, int, int]] = None) -> None:
        r2 = radius * radius
        r_in2 = max(0, (radius - 1) * (radius - 1))
        for dy in range(-radius, radius + 1):
            for dx in range(-radius, radius + 1):
                d2 = dx * dx + dy * dy
                if d2 <= r2:
                    if outline and d2 >= r_in2:
                        set_px(px + dx, py + dy, *outline)
                    else:
                        set_px(px + dx, py + dy, *rgb)

    def draw_line(x0: int, y0: int, x1: int, y1: int, rgb: Tuple[int, int, int], thickness: int = 2) -> None:
        steps = max(abs(x1 - x0), abs(y1 - y0), 1)
        for i in range(steps + 1):
            x = int(round(x0 + (x1 - x0) * i / steps))
            y = int(round(y0 + (y1 - y0) * i / steps))
            draw_circle(x, y, thickness, rgb)

    # スケール計算
    all_pts = [
        rot(a.x - cx, a.y - cy, a.z - cz)
        for residues in protein_chains.values()
        for res in residues
        for a in res.heavy_atoms
    ]
    max_span = max([max(abs(p[0]), abs(p[1])) for p in all_pts] + [10.0])
    scale = min(width, height) * 0.38 / max_span

    def to_screen(x: float, y: float, z: float) -> Tuple[int, int]:
        rx, ry, _ = rot(x - cx, y - cy, z - cz)
        return (int(width / 2 + rx * scale), int(height / 2 - ry * scale))

    # 外枠
    for x in range(width):
        set_px(x, 0, 200, 200, 200)
        set_px(x, height - 1, 200, 200, 200)
    for y in range(height):
        set_px(0, y, 200, 200, 200)
        set_px(width - 1, y, 200, 200, 200)

    # 1. 酵素主鎖 (青色: #1565C0)
    for residues in protein_chains.values():
        pts = [to_screen(r.ca_atom.x, r.ca_atom.y, r.ca_atom.z) for r in residues]
        for i in range(len(pts) - 1):
            draw_line(pts[i][0], pts[i][1], pts[i + 1][0], pts[i + 1][1], (25, 100, 210), thickness=2)
        for px, py in pts:
            draw_circle(px, py, 3, (35, 115, 225))

    # 2. リガンド: 基質・補酵素 (黄色: #FFD600)
    for lig in ligands:
        l_pts = [to_screen(a.x, a.y, a.z) for a in lig.heavy_atoms]
        for i in range(len(l_pts) - 1):
            draw_line(l_pts[i][0], l_pts[i][1], l_pts[i + 1][0], l_pts[i + 1][1], (230, 180, 0), thickness=2)
        for px, py in l_pts:
            draw_circle(px, py, 5, (255, 215, 0), outline=(140, 105, 0))

    # 3. N末端 (赤色: #E50000) & C末端 (緑色: #00B840)
    for residues in protein_chains.values():
        n_pt = to_screen(residues[0].ca_atom.x, residues[0].ca_atom.y, residues[0].ca_atom.z)
        c_pt = to_screen(residues[-1].ca_atom.x, residues[-1].ca_atom.y, residues[-1].ca_atom.z)
        draw_circle(n_pt[0], n_pt[1], 9, (230, 20, 20), outline=(40, 0, 0))
        draw_circle(c_pt[0], c_pt[1], 9, (0, 185, 65), outline=(0, 40, 10))

    # PNG エンコード
    raw_rows = bytearray()
    stride = width * 3
    for y in range(height):
        raw_rows.append(0)  # Filter type 0 (None)
        raw_rows.extend(buf[y * stride : (y + 1) * stride])

    def png_chunk(chunk_type: bytes, data: bytes) -> bytes:
        crc = zlib.crc32(chunk_type + data) & 0xFFFFFFFF
        return struct.pack(">I", len(data)) + chunk_type + data + struct.pack(">I", crc)

    ihdr = struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)
    png_bytes = (
        b"\x89PNG\r\n\x1a\n"
        + png_chunk(b"IHDR", ihdr)
        + png_chunk(b"IDAT", zlib.compress(bytes(raw_rows), 9))
        + png_chunk(b"IEND", b"")
    )
    out_path.write_bytes(png_bytes)


def render_three_views_for_pdb(
    result: PDBEvaluationResult,
    output_img_dir: Path,
    pymol_exe: Optional[str],
    rotation_mode: str = "y90_x90",
    image_width: int = 1200,
    image_height: int = 900,
    ray_trace: bool = True,
    allow_fallback: bool = True,
) -> List[Path]:
    """
    1つのPDB構造に対して、PyMOL を用いて3方向（0°、90°回転#1、90°回転#2）のPNG画像を作成し、
    画像パスのリストを返す。また再現用 .pml スクリプトも保存する。
    """
    output_img_dir.mkdir(parents=True, exist_ok=True)
    pml_dir = output_img_dir / "pml_scripts"
    pml_dir.mkdir(parents=True, exist_ok=True)

    stem = result.pdb_name
    img1 = output_img_dir / f"{stem}_view1_0deg.png"
    img2 = output_img_dir / f"{stem}_view2_rot90_1.png"
    img3 = output_img_dir / f"{stem}_view3_rot90_2.png"
    out_paths = [img1, img2, img3]

    # WSL2 から Windows の PyMOL.exe を呼ぶ場合、パスが /home/... だと Windows 側から直接開けない場合があるため、
    # 必要に応じて Windows 側の一時ディレクトリまたは変換パスを使用する
    is_win_exe_on_wsl = bool(is_wsl2() and pymol_exe and pymol_exe.lower().endswith(".exe"))

    pdb_for_pml = to_pymol_compatible_path(result.pdb_path, pymol_exe)
    img1_for_pml = to_pymol_compatible_path(img1, pymol_exe)
    img2_for_pml = to_pymol_compatible_path(img2, pymol_exe)
    img3_for_pml = to_pymol_compatible_path(img3, pymol_exe)

    pml_content = build_pymol_script_content(
        pdb_file_for_pymol=pdb_for_pml,
        out_img1=img1_for_pml,
        out_img2=img2_for_pml,
        out_img3=img3_for_pml,
        result=result,
        rotation_mode=rotation_mode,
        image_width=image_width,
        image_height=image_height,
        ray_trace=ray_trace,
    )

    pml_path = pml_dir / f"{stem}_render.pml"
    pml_path.write_text(pml_content, encoding="utf-8")
    result.pml_path = pml_path

    if pymol_exe:
        try:
            if is_win_exe_on_wsl and not str(result.pdb_path.resolve()).startswith("/mnt/"):
                # WSL Linux ファイルシステム (/home/...) にある場合、/mnt/c/Users/Public または Windows Temp にステージング
                win_temp_candidates = glob.glob("/mnt/c/Users/*/AppData/Local/Temp")
                stage_root = Path(win_temp_candidates[0]) if win_temp_candidates else Path("/mnt/c/Temp")
                stage_dir = stage_root / f"pymol_histag_{os.getpid()}"
                stage_dir.mkdir(parents=True, exist_ok=True)
                try:
                    staged_pdb = stage_dir / result.pdb_path.name
                    shutil.copy2(result.pdb_path, staged_pdb)
                    s_img1 = stage_dir / img1.name
                    s_img2 = stage_dir / img2.name
                    s_img3 = stage_dir / img3.name
                    staged_pml = stage_dir / pml_path.name
                    staged_pml.write_text(
                        build_pymol_script_content(
                            pdb_file_for_pymol=to_pymol_compatible_path(staged_pdb, pymol_exe),
                            out_img1=to_pymol_compatible_path(s_img1, pymol_exe),
                            out_img2=to_pymol_compatible_path(s_img2, pymol_exe),
                            out_img3=to_pymol_compatible_path(s_img3, pymol_exe),
                            result=result,
                            rotation_mode=rotation_mode,
                            image_width=image_width,
                            image_height=image_height,
                            ray_trace=ray_trace,
                        ),
                        encoding="utf-8",
                    )
                    pml_arg = to_pymol_compatible_path(staged_pml, pymol_exe)
                    subprocess.run([pymol_exe, "-cq", pml_arg], check=True, capture_output=True, text=True)
                    for src, dst in ((s_img1, img1), (s_img2, img2), (s_img3, img3)):
                        if src.is_file():
                            shutil.copy2(src, dst)
                finally:
                    shutil.rmtree(stage_dir, ignore_errors=True)
            else:
                pml_arg = to_pymol_compatible_path(pml_path, pymol_exe)
                subprocess.run([pymol_exe, "-cq", pml_arg], check=True, capture_output=True, text=True)

            if all(p.is_file() for p in out_paths):
                result.image_paths = out_paths
                return out_paths
        except Exception as exc:
            print(f"[WARNING] PyMOL 実行中に警告が発生しました ({stem}): {exc}", file=sys.stderr)

    if allow_fallback:
        _render_fallback_3d_images(result.pdb_path, result, out_paths, rotation_mode=rotation_mode)
        result.image_paths = out_paths

    return out_paths


# =============================================================================
# 7. Excel (.xlsx) レポート出力 (画像3種貼り付け + 論文基準数値 + 参考文献シート)
# =============================================================================

REFERENCES_DATA: List[Tuple[str, str, str, str]] = [
    (
        "[1] SASA (溶媒露出表面積) 計算アルゴリズム",
        "Shrake, A., & Rupley, J. A. (1973). Environment and exposure to solvent of protein atoms. Lysozyme and insulin. Journal of Molecular Biology, 79(2), 351-371.",
        "10.1016/0022-2836(73)90011-9",
        "水分子プローブ半径 1.4 Å を用いた Rolling-ball (球面点群) 法により各原子・残基の絶対 SASA (Å²) を算出。",
    ),
    (
        "[2] RSA (相対溶媒露出度) 正規化基準値 (MaxASA)",
        "Tien, M. Z., Meyer, A. G., Sydykova, D. K., Spielman, S. J., & Wilke, C. O. (2013). Maximum allowed solvent accessibilities of residues in proteins. PLOS ONE, 8(11), e80635.",
        "10.1371/journal.pone.0080635",
        "Gly-X-Gly トリペプチドの理論的最大露出表面積 (MaxASA) により SASA を正規化し、RSA = SASA / MaxASA (0〜1) を算出。",
    ),
    (
        "[3] タンパク質 N末端・C末端の溶媒露出性と構造的役割",
        "Thornton, J. M., & Sibanda, B. L. (1983). Amino and carboxy-terminal regions in globular proteins. Journal of Molecular Biology, 167(2), 443-460.",
        "10.1016/S0022-2836(83)80344-1",
        "球状タンパク質における N末・C末領域の露出度・柔軟性・空間配置の統計解析。最末端1残基だけでなく末端近傍窓 (N1–N3, C1–C3) の平均露出度評価が有効であることを示す。",
    ),
    (
        "[4] 埋没残基・露出残基の RSA 閾値分類 (1)",
        "Miller, S., Janin, J., Lesk, A. M., & Chothia, C. (1987). Interior and surface of monomeric proteins. Journal of Molecular Biology, 196(3), 641-656.",
        "10.1016/0022-2836(87)90038-6",
        "RSA < 0.20 (20%) を疎水性コアへの埋没残基 (Buried) と定義する標準的クライテリア。",
    ),
    (
        "[5] 埋没残基・露出残基の RSA 閾値分類 (2)",
        "Rost, B., & Sander, C. (1994). Conservation and prediction of solvent accessibility in protein families. Proteins: Structure, Function, and Bioinformatics, 20(3), 216-226.",
        "10.1002/prot.340200303",
        "RSA < 0.20 (埋没: 高リスク), 0.20 <= RSA < 0.50 (部分露出: 中リスク), RSA >= 0.50 (十分露出: 低リスク) の3段階評価基準。",
    ),
    (
        "[6] Hisタグの立体構造への影響と物理的リーチ長",
        "Carson, M., Johnson, D. H., McDonald, H., Brouillette, C., & DeLucas, L. J. (2007). His-tag impact on structure. Acta Crystallographica Section D, 63(3), 295-301.",
        "10.1107/S0907444906052024",
        "PDB内のHisタグ付加構造の網羅的解析。ヘキサヒスチジン (His6) の完全伸長時の物理的コンター長は約 3.6 Å/残基 × 6 ≈ 21.6 Å（ランダムコイル実効到達半径 ~12–15 Å）。",
    ),
    (
        "[7] His6タグと酵素活性部位の立体的・静電的干渉モデリング",
        "Freydank, A. C., Brandt, W., & Dräger, B. (2008). Protein structure modeling indicates hexahistidine-tag interference with enzyme activity. Proteins, 72(1), 173-183.",
        "10.1002/prot.21905",
        "末端に付加した His6 タグが活性部位ポケットに到達可能な距離にある場合、基質結合の立体障害や活性低下を引き起こすことを構造モデリングにより実証。",
    ),
    (
        "[8] N末端 vs C末端 Hisタグ配置による酵素触媒活性への影響",
        "Yeon, Y. J., Park, H. J., & Yoo, Y. J. (2014). Effect of His-tag location on the catalytic activity of 3-hydroxybutyrate dehydrogenase. Biotechnology and Bioprocess Engineering, 19(5), 798-802.",
        "10.1007/s12257-014-0089-2",
        "基質結合ポケット・ヒンジ領域に近い側の末端へのHisタグ付加が酵素の触媒効率を著しく低下させる一方、活性部位から離れた末端への付加では活性が維持されることを実証。",
    ),
]


def export_to_excel(results: List[PDBEvaluationResult], output_excel_path: Path) -> None:
    """
    評価結果と3方向の構造図を Excel (.xlsx) ファイルに出力する。
    openpyxl がインストールされていれば openpyxl を使用し、未インストールの場合は
    標準ライブラリのみで画像を埋め込んだ正規の OpenXML (.xlsx) ファイルを出力する。
    """
    output_excel_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        import openpyxl  # noqa: F401
        _export_to_excel_openpyxl(results, output_excel_path)
    except ImportError:
        _export_to_excel_stdlib(results, output_excel_path)


def _export_to_excel_openpyxl(results: List[PDBEvaluationResult], output_excel_path: Path) -> None:
    import openpyxl
    from openpyxl.drawing.image import Image as XLImage
    from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
    from openpyxl.utils import get_column_letter

    wb = openpyxl.Workbook()

    # --- Sheet 1: Summary & 3-View Structures ---
    ws1 = wb.active
    ws1.title = "Summary_&_3Views"

    header_fill = PatternFill("solid", fgColor="1F4E78")
    header_font = Font(color="FFFFFF", bold=True, size=10)
    thin_side = Side(border_style="thin", color="CCCCCC")
    border = Border(left=thin_side, right=thin_side, top=thin_side, bottom=thin_side)

    fill_high = PatternFill("solid", fgColor="FCE4D6")
    fill_mod = PatternFill("solid", fgColor="FFF2CC")
    fill_low = PatternFill("solid", fgColor="E2EFDA")

    headers_s1 = [
        "PDB Name",
        "Chains / Residues",
        "Detected Ligands\n(Substrate & Cofactor)",
        "Recommended\nHis-Tag Position",
        # N-terminus
        "N-Term\nResidue",
        "N-Term Min Dist\nto All Ligands (Å)",
        "N-Term Dist to\nSubstrate (Å)",
        "N-Term Dist to\nCofactor (Å)",
        "N-Term SASA\n(1-res / 3-win, Å²)",
        "N-Term RSA\n(1-res / 3-win)",
        "N-Term Exposure\nClass",
        "N-Term Risk\nScore (0-100)",
        "N-Term Risk\nLevel",
        # C-terminus
        "C-Term\nResidue",
        "C-Term Min Dist\nto All Ligands (Å)",
        "C-Term Dist to\nSubstrate (Å)",
        "C-Term Dist to\nCofactor (Å)",
        "C-Term SASA\n(1-res / 3-win, Å²)",
        "C-Term RSA\n(1-res / 3-win)",
        "C-Term Exposure\nClass",
        "C-Term Risk\nScore (0-100)",
        "C-Term Risk\nLevel",
        # Design Note & 3 Images
        "Design Assessment &\nLinker Recommendation",
        "View 1 (0° Oriented)\n[Blue:Enzyme, Yellow:Ligands,\nRed:N-term, Green:C-term]",
        "View 2 (+90° Rot #1)\n[Blue:Enzyme, Yellow:Ligands,\nRed:N-term, Green:C-term]",
        "View 3 (+90° Rot #2)\n[Blue:Enzyme, Yellow:Ligands,\nRed:N-term, Green:C-term]",
    ]

    ws1.append(headers_s1)
    ws1.row_dimensions[1].height = 42
    for col_idx in range(1, len(headers_s1) + 1):
        cell = ws1.cell(row=1, column=col_idx)
        cell.fill = header_fill
        cell.font = header_font
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        cell.border = border

    img_col_start = 24  # Columns X, Y, Z (24, 25, 26)
    for row_idx, res in enumerate(results, start=2):
        lig_summary = ", ".join([f"{l.role}:{l.label}" for l in res.ligands]) if res.ligands else "None"
        nt = res.n_term
        ct = res.c_term

        row_vals = [
            res.pdb_name,
            f"{','.join(res.protein_chains)} ({res.total_residues} aa)",
            lig_summary,
            res.recommended_terminus,
            # N-term
            f"{nt.chain}:{nt.residue_label}",
            nt.dist_min_all_ligands if nt.dist_min_all_ligands is not None else "N/A",
            nt.dist_min_substrate if nt.dist_min_substrate is not None else "N/A",
            nt.dist_min_cofactor if nt.dist_min_cofactor is not None else "N/A",
            f"{nt.sasa_1res:.1f} / {nt.sasa_win3:.1f}",
            f"{nt.rsa_1res*100:.1f}% / {nt.rsa_win3*100:.1f}%",
            nt.exposure_class,
            nt.total_risk_score,
            nt.risk_level,
            # C-term
            f"{ct.chain}:{ct.residue_label}",
            ct.dist_min_all_ligands if ct.dist_min_all_ligands is not None else "N/A",
            ct.dist_min_substrate if ct.dist_min_substrate is not None else "N/A",
            ct.dist_min_cofactor if ct.dist_min_cofactor is not None else "N/A",
            f"{ct.sasa_1res:.1f} / {ct.sasa_win3:.1f}",
            f"{ct.rsa_1res*100:.1f}% / {ct.rsa_win3*100:.1f}%",
            ct.exposure_class,
            ct.total_risk_score,
            ct.risk_level,
            # Comment + placeholders for images
            f"{res.recommendation_reason}\n[N末] {nt.design_comment}\n[C末] {ct.design_comment}",
            "",
            "",
            "",
        ]
        ws1.append(row_vals)
        ws1.row_dimensions[row_idx].height = 145

        for c_i in range(1, len(headers_s1) + 1):
            c = ws1.cell(row=row_idx, column=c_i)
            c.alignment = Alignment(
                horizontal="left" if c_i in (3, 23) else "center",
                vertical="center",
                wrap_text=True,
            )
            c.border = border

        # リスクレベルの色分け
        for col_risk, level_str in ((12, nt.risk_level), (13, nt.risk_level), (21, ct.risk_level), (22, ct.risk_level)):
            cell = ws1.cell(row=row_idx, column=col_risk)
            if level_str.startswith("HIGH"):
                cell.fill = fill_high
            elif level_str.startswith("MODERATE"):
                cell.fill = fill_mod
            else:
                cell.fill = fill_low

        # 3方向の画像をセルに埋め込み
        for v_i, img_path in enumerate(res.image_paths[:3]):
            if img_path and img_path.is_file():
                xl_img = XLImage(str(img_path))
                xl_img.width = 230
                xl_img.height = 172
                col_letter = get_column_letter(img_col_start + v_i)
                ws1.add_image(xl_img, f"{col_letter}{row_idx}")

    col_widths = {
        1: 18, 2: 14, 3: 24, 4: 24,
        5: 12, 6: 15, 7: 15, 8: 15, 9: 16, 10: 16, 11: 22, 12: 13, 13: 18,
        14: 12, 15: 15, 16: 15, 17: 15, 18: 16, 19: 16, 20: 22, 21: 13, 22: 18,
        23: 42, 24: 34, 25: 34, 26: 34,
    }
    for c_i, w in col_widths.items():
        ws1.column_dimensions[get_column_letter(c_i)].width = w

    # --- Sheet 2: Detailed Terminus Metrics ---
    ws2 = wb.create_sheet(title="Detailed_Metrics")
    headers_s2 = [
        "PDB Name", "Chain", "Terminus", "Residue",
        "SASA 1-res (Å²)", "RSA 1-res (Tien 2013)",
        "SASA 3-res Window (Å²)", "RSA 3-res Window (Thornton 1983 / Tien 2013)",
        "Exposure Classification", "Mean B-factor / pLDDT",
        "Nearest Ligand", "Min Heavy-Atom Dist (All Ligands, Å)", "Cα-Centroid Dist (All Ligands, Å)",
        "Substrate ID", "Min Dist to Substrate (Å)", "Cα-Centroid Dist to Substrate (Å)",
        "Cofactor ID", "Min Dist to Cofactor (Å)", "Cα-Centroid Dist to Cofactor (Å)",
        "Terminal Orientation Angle θ (°)",
        "Distance Risk Score", "RSA Burial Risk Score", "Orientation Risk Score",
        "Total Risk Score (0-100)", "Risk Level", "Design Notes",
    ]
    ws2.append(headers_s2)
    ws2.row_dimensions[1].height = 32
    for c_i in range(1, len(headers_s2) + 1):
        c = ws2.cell(row=1, column=c_i)
        c.fill = header_fill
        c.font = header_font
        c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        c.border = border

    r2 = 2
    for res in results:
        for n_m, c_m in res.all_chain_termini:
            for tm in (n_m, c_m):
                ws2.append([
                    res.pdb_name, tm.chain, tm.terminus_type, tm.residue_label,
                    tm.sasa_1res, tm.rsa_1res, tm.sasa_win3, tm.rsa_win3,
                    tm.exposure_class, tm.mean_bfactor,
                    tm.nearest_ligand_label,
                    tm.dist_min_all_ligands if tm.dist_min_all_ligands is not None else "N/A",
                    tm.dist_ca_centroid_all if tm.dist_ca_centroid_all is not None else "N/A",
                    tm.substrate_label,
                    tm.dist_min_substrate if tm.dist_min_substrate is not None else "N/A",
                    tm.dist_ca_substrate if tm.dist_ca_substrate is not None else "N/A",
                    tm.cofactor_label,
                    tm.dist_min_cofactor if tm.dist_min_cofactor is not None else "N/A",
                    tm.dist_ca_cofactor if tm.dist_ca_cofactor is not None else "N/A",
                    tm.orientation_angle_deg if tm.orientation_angle_deg is not None else "N/A",
                    tm.dist_risk_score, tm.rsa_risk_score, tm.orient_risk_score,
                    tm.total_risk_score, tm.risk_level, tm.design_comment,
                ])
                for c_i in range(1, len(headers_s2) + 1):
                    cell = ws2.cell(row=r2, column=c_i)
                    cell.border = border
                    cell.alignment = Alignment(horizontal="left" if c_i == 26 else "center", vertical="center")
                r2 += 1

    for c_i in range(1, len(headers_s2) + 1):
        ws2.column_dimensions[get_column_letter(c_i)].width = 38 if c_i == 26 else 18

    # --- Sheet 3: Methodology & Literature References ---
    ws3 = wb.create_sheet(title="Methodology_&_References")
    ws3.append(["評価項目 / 根拠カテゴリ", "引用論文 (Citation)", "DOI", "本スクリプトにおける評価基準・数式への適用内容"])
    ws3.row_dimensions[1].height = 28
    for c_i in range(1, 5):
        c = ws3.cell(row=1, column=c_i)
        c.fill = header_fill
        c.font = header_font
        c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        c.border = border

    for idx, (cat, cit, doi, desc) in enumerate(REFERENCES_DATA, start=2):
        ws3.append([cat, cit, doi, desc])
        ws3.row_dimensions[idx].height = 36
        for c_i in range(1, 5):
            c = ws3.cell(row=idx, column=c_i)
            c.alignment = Alignment(horizontal="left", vertical="center", wrap_text=True)
            c.border = border

    ws3.column_dimensions["A"].width = 34
    ws3.column_dimensions["B"].width = 65
    ws3.column_dimensions["C"].width = 28
    ws3.column_dimensions["D"].width = 65

    wb.save(str(output_excel_path))


def _export_to_excel_stdlib(results: List[PDBEvaluationResult], output_excel_path: Path) -> None:
    """
    openpyxl が未インストールの環境でも動作する、標準ライブラリ (zipfile + XML) による
    画像付き .xlsx (Office Open XML) エクスポーター。
    3つのシート (Summary_&_3Views, Detailed_Metrics, Methodology_&_References) と
    各PDBにつき3枚のPNG画像を埋め込んだ正規の Excel ワークブックを生成する。
    """
    # Sheet 1 の行データを構築
    headers_s1 = [
        "PDB Name", "Chains / Residues", "Detected Ligands (Substrate & Cofactor)",
        "Recommended His-Tag Position",
        "N-Term Residue", "N-Term Min Dist All Ligands (Å)", "N-Term Dist Substrate (Å)",
        "N-Term Dist Cofactor (Å)", "N-Term SASA 1res/3win (Å²)", "N-Term RSA 1res/3win",
        "N-Term Exposure Class", "N-Term Risk Score (0-100)", "N-Term Risk Level",
        "C-Term Residue", "C-Term Min Dist All Ligands (Å)", "C-Term Dist Substrate (Å)",
        "C-Term Dist Cofactor (Å)", "C-Term SASA 1res/3win (Å²)", "C-Term RSA 1res/3win",
        "C-Term Exposure Class", "C-Term Risk Score (0-100)", "C-Term Risk Level",
        "Design Assessment & Linker Recommendation",
        "View 1 (0 deg Oriented)", "View 2 (+90 deg Rot #1)", "View 3 (+90 deg Rot #2)",
    ]
    rows_s1: List[List[object]] = [headers_s1]
    media_files: List[Tuple[str, Path, int, int]] = []  # (media_name, path, row_0based, col_0based)

    for r_idx, res in enumerate(results, start=1):
        lig_summary = ", ".join([f"{l.role}:{l.label}" for l in res.ligands]) if res.ligands else "None"
        nt, ct = res.n_term, res.c_term
        rows_s1.append([
            res.pdb_name,
            f"{','.join(res.protein_chains)} ({res.total_residues} aa)",
            lig_summary,
            res.recommended_terminus,
            f"{nt.chain}:{nt.residue_label}",
            nt.dist_min_all_ligands if nt.dist_min_all_ligands is not None else "N/A",
            nt.dist_min_substrate if nt.dist_min_substrate is not None else "N/A",
            nt.dist_min_cofactor if nt.dist_min_cofactor is not None else "N/A",
            f"{nt.sasa_1res:.1f} / {nt.sasa_win3:.1f}",
            f"{nt.rsa_1res*100:.1f}% / {nt.rsa_win3*100:.1f}%",
            nt.exposure_class,
            nt.total_risk_score,
            nt.risk_level,
            f"{ct.chain}:{ct.residue_label}",
            ct.dist_min_all_ligands if ct.dist_min_all_ligands is not None else "N/A",
            ct.dist_min_substrate if ct.dist_min_substrate is not None else "N/A",
            ct.dist_min_cofactor if ct.dist_min_cofactor is not None else "N/A",
            f"{ct.sasa_1res:.1f} / {ct.sasa_win3:.1f}",
            f"{ct.rsa_1res*100:.1f}% / {ct.rsa_win3*100:.1f}%",
            ct.exposure_class,
            ct.total_risk_score,
            ct.risk_level,
            f"{res.recommendation_reason} | [N末] {nt.design_comment} | [C末] {ct.design_comment}",
            "", "", "",
        ])
        for v_i, img_p in enumerate(res.image_paths[:3]):
            if img_p and img_p.is_file():
                media_name = f"image_{r_idx}_{v_i + 1}.png"
                media_files.append((media_name, img_p, r_idx, 23 + v_i))

    # Sheet 2
    headers_s2 = [
        "PDB Name", "Chain", "Terminus", "Residue",
        "SASA 1-res (Å²)", "RSA 1-res (Tien 2013)",
        "SASA 3-res Window (Å²)", "RSA 3-res Window (Thornton 1983 / Tien 2013)",
        "Exposure Classification", "Mean B-factor / pLDDT",
        "Nearest Ligand", "Min Heavy-Atom Dist (All Ligands, Å)", "Cα-Centroid Dist (All Ligands, Å)",
        "Substrate ID", "Min Dist to Substrate (Å)", "Cα-Centroid Dist to Substrate (Å)",
        "Cofactor ID", "Min Dist to Cofactor (Å)", "Cα-Centroid Dist to Cofactor (Å)",
        "Terminal Orientation Angle θ (°)",
        "Distance Risk Score", "RSA Burial Risk Score", "Orientation Risk Score",
        "Total Risk Score (0-100)", "Risk Level", "Design Notes",
    ]
    rows_s2: List[List[object]] = [headers_s2]
    for res in results:
        for n_m, c_m in res.all_chain_termini:
            for tm in (n_m, c_m):
                rows_s2.append([
                    res.pdb_name, tm.chain, tm.terminus_type, tm.residue_label,
                    tm.sasa_1res, tm.rsa_1res, tm.sasa_win3, tm.rsa_win3,
                    tm.exposure_class, tm.mean_bfactor,
                    tm.nearest_ligand_label,
                    tm.dist_min_all_ligands if tm.dist_min_all_ligands is not None else "N/A",
                    tm.dist_ca_centroid_all if tm.dist_ca_centroid_all is not None else "N/A",
                    tm.substrate_label,
                    tm.dist_min_substrate if tm.dist_min_substrate is not None else "N/A",
                    tm.dist_ca_substrate if tm.dist_ca_substrate is not None else "N/A",
                    tm.cofactor_label,
                    tm.dist_min_cofactor if tm.dist_min_cofactor is not None else "N/A",
                    tm.dist_ca_cofactor if tm.dist_ca_cofactor is not None else "N/A",
                    tm.orientation_angle_deg if tm.orientation_angle_deg is not None else "N/A",
                    tm.dist_risk_score, tm.rsa_risk_score, tm.orient_risk_score,
                    tm.total_risk_score, tm.risk_level, tm.design_comment,
                ])

    # Sheet 3
    rows_s3: List[List[object]] = [
        ["評価項目 / 根拠カテゴリ", "引用論文 (Citation)", "DOI", "本スクリプトにおける評価基準・数式への適用内容"]
    ]
    for cat, cit, doi, desc in REFERENCES_DATA:
        rows_s3.append([cat, cit, doi, desc])

    def col_name(col_idx_0: int) -> str:
        name = ""
        n = col_idx_0 + 1
        while n > 0:
            n, rem = divmod(n - 1, 26)
            name = chr(65 + rem) + name
        return name

    def build_sheet_xml(rows: List[List[object]], row_height_data: float = 24.0, has_drawing: bool = False) -> str:
        lines = [
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
            '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
            'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">',
            "<sheetData>",
        ]
        for r_i, row in enumerate(rows, start=1):
            ht = 32.0 if r_i == 1 else row_height_data
            lines.append(f'<row r="{r_i}" ht="{ht}" customHeight="1">')
            for c_i, val in enumerate(row):
                cell_ref = f"{col_name(c_i)}{r_i}"
                style_attr = ' s="1"' if r_i == 1 else ' s="2"'
                if isinstance(val, (int, float)) and not isinstance(val, bool):
                    lines.append(f'<c r="{cell_ref}"{style_attr}><v>{val}</v></c>')
                else:
                    s_val = xml_escape(str(val))
                    lines.append(f'<c r="{cell_ref}" t="inlineStr"{style_attr}><is><t>{s_val}</t></is></c>')
            lines.append("</row>")
        lines.append("</sheetData>")
        if has_drawing:
            lines.append('<drawing r:id="rId1"/>')
        lines.append("</worksheet>")
        return "\n".join(lines)

    # Drawing XML (画像のセルアンカー定義)
    drawing_anchors: List[str] = []
    drawing_rels: List[str] = []
    for idx, (media_name, _, r_0, c_0) in enumerate(media_files, start=1):
        rid = f"rId{idx}"
        drawing_Anch = f"""
  <xdr:twoCellAnchor editAs="oneCell">
    <xdr:from><xdr:col>{c_0}</xdr:col><xdr:colOff>19050</xdr:colOff><xdr:row>{r_0}</xdr:row><xdr:rowOff>19050</xdr:rowOff></xdr:from>
    <xdr:to><xdr:col>{c_0 + 1}</xdr:col><xdr:colOff>-19050</xdr:colOff><xdr:row>{r_0 + 1}</xdr:row><xdr:rowOff>-19050</xdr:rowOff></xdr:to>
    <xdr:pic>
      <xdr:nvPicPr>
        <xdr:cNvPr id="{idx}" name="Picture {idx}"/>
        <xdr:cNvPicPr><a:picLocks noChangeAspect="1"/></xdr:cNvPicPr>
      </xdr:nvPicPr>
      <xdr:blipFill>
        <a:blip r:embed="{rid}"/>
        <a:stretch><a:fillRect/></a:stretch>
      </xdr:blipFill>
      <xdr:spPr>
        <a:xfrm><a:off x="0" y="0"/><a:ext cx="2190750" cy="1638300"/></a:xfrm>
        <a:prstGeom prst="rect"><a:avLst/></a:prstGeom>
      </xdr:spPr>
    </xdr:pic>
    <xdr:clientData/>
  </xdr:twoCellAnchor>"""
        drawing_anchors.append(drawing_Anch)
        drawing_rels.append(
            f'<Relationship Id="{rid}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="../media/{media_name}"/>'
        )

    drawing_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
        '<xdr:wsDr xmlns:xdr="http://schemas.openxmlformats.org/drawingml/2006/spreadsheetDrawing" '
        'xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">\n'
        + "\n".join(drawing_anchors)
        + "\n</xdr:wsDr>"
    )
    drawing_rels_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">\n'
        + "\n".join(drawing_rels)
        + "\n</Relationships>"
    )

    styles_xml = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
  <fonts count="2">
    <font><sz val="10"/><name val="Calibri"/></font>
    <font><b/><sz val="10"/><color rgb="FFFFFFFF"/><name val="Calibri"/></font>
  </fonts>
  <fills count="3">
    <fill><patternFill patternType="none"/></fill>
    <fill><patternFill patternType="gray125"/></fill>
    <fill><patternFill patternType="solid"><fgColor rgb="FF1F4E78"/><bgColor indexed="64"/></patternFill></fill>
  </fills>
  <borders count="2">
    <border><left/><right/><top/><bottom/><diagonal/></border>
    <border>
      <left style="thin"><color rgb="FFCCCCCC"/></left>
      <right style="thin"><color rgb="FFCCCCCC"/></right>
      <top style="thin"><color rgb="FFCCCCCC"/></top>
      <bottom style="thin"><color rgb="FFCCCCCC"/></bottom>
      <diagonal/>
    </border>
  </borders>
  <cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>
  <cellXfs count="3">
    <xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>
    <xf numFmtId="0" fontId="1" fillId="2" borderId="1" xfId="0" applyFont="1" applyFill="1" applyBorder="1" applyAlignment="1">
      <alignment horizontal="center" vertical="center" wrapText="1"/>
    </xf>
    <xf numFmtId="0" fontId="0" fillId="0" borderId="1" xfId="0" applyBorder="1" applyAlignment="1">
      <alignment vertical="center" wrapText="1"/>
    </xf>
  </cellXfs>
</styleSheet>"""

    workbook_xml = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
          xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <sheets>
    <sheet name="Summary_&amp;_3Views" sheetId="1" r:id="rId1"/>
    <sheet name="Detailed_Metrics" sheetId="2" r:id="rId2"/>
    <sheet name="Methodology_&amp;_References" sheetId="3" r:id="rId3"/>
  </sheets>
</workbook>"""

    workbook_rels_xml = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet2.xml"/>
  <Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet3.xml"/>
  <Relationship Id="rId4" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>
</Relationships>"""

    root_rels_xml = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
</Relationships>"""

    sheet1_rels_xml = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/drawing" Target="../drawings/drawing1.xml"/>
</Relationships>"""

    content_types_xml = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Default Extension="png" ContentType="image/png"/>
  <Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
  <Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>
  <Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>
  <Override PartName="/xl/worksheets/sheet2.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>
  <Override PartName="/xl/worksheets/sheet3.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>
  <Override PartName="/xl/drawings/drawing1.xml" ContentType="application/vnd.openxmlformats-officedocument.drawing+xml"/>
</Types>"""

    has_imgs = len(media_files) > 0
    with zipfile.ZipFile(output_excel_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", content_types_xml)
        zf.writestr("_rels/.rels", root_rels_xml)
        zf.writestr("xl/workbook.xml", workbook_xml)
        zf.writestr("xl/_rels/workbook.xml.rels", workbook_rels_xml)
        zf.writestr("xl/styles.xml", styles_xml)
        zf.writestr("xl/worksheets/sheet1.xml", build_sheet_xml(rows_s1, row_height_data=145.0, has_drawing=has_imgs))
        zf.writestr("xl/worksheets/sheet2.xml", build_sheet_xml(rows_s2, row_height_data=24.0, has_drawing=False))
        zf.writestr("xl/worksheets/sheet3.xml", build_sheet_xml(rows_s3, row_height_data=34.0, has_drawing=False))
        if has_imgs:
            zf.writestr("xl/worksheets/_rels/sheet1.xml.rels", sheet1_rels_xml)
            zf.writestr("xl/drawings/drawing1.xml", drawing_xml)
            zf.writestr("xl/drawings/_rels/drawing1.xml.rels", drawing_rels_xml)
            for media_name, img_p, _, _ in media_files:
                zf.write(img_p, f"xl/media/{media_name}")


def export_to_csv(results: List[PDBEvaluationResult], output_csv_path: Path) -> None:
    """Excelファイルと併せて、機械可読なサマリー CSV ファイルも出力する。"""
    output_csv_path.parent.mkdir(parents=True, exist_ok=True)
    headers = [
        "pdb_name", "chains", "total_residues", "ligands", "recommended_terminus",
        "n_term_res", "n_term_dist_min_all_A", "n_term_dist_substrate_A", "n_term_dist_cofactor_A",
        "n_term_sasa_1res_A2", "n_term_rsa_1res", "n_term_sasa_win3_A2", "n_term_rsa_win3",
        "n_term_exposure_class", "n_term_orientation_deg", "n_term_risk_score", "n_term_risk_level",
        "c_term_res", "c_term_dist_min_all_A", "c_term_dist_substrate_A", "c_term_dist_cofactor_A",
        "c_term_sasa_1res_A2", "c_term_rsa_1res", "c_term_sasa_win3_A2", "c_term_rsa_win3",
        "c_term_exposure_class", "c_term_orientation_deg", "c_term_risk_score", "c_term_risk_level",
        "view1_image", "view2_image", "view3_image",
    ]
    with open(output_csv_path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        for res in results:
            nt, ct = res.n_term, res.c_term
            imgs = [str(p) for p in res.image_paths] + ["", "", ""]
            writer.writerow([
                res.pdb_name,
                ";".join(res.protein_chains),
                res.total_residues,
                ";".join([f"{l.role}:{l.label}" for l in res.ligands]),
                res.recommended_terminus,
                f"{nt.chain}:{nt.residue_label}",
                nt.dist_min_all_ligands, nt.dist_min_substrate, nt.dist_min_cofactor,
                nt.sasa_1res, nt.rsa_1res, nt.sasa_win3, nt.rsa_win3,
                nt.exposure_class, nt.orientation_angle_deg, nt.total_risk_score, nt.risk_level,
                f"{ct.chain}:{ct.residue_label}",
                ct.dist_min_all_ligands, ct.dist_min_substrate, ct.dist_min_cofactor,
                ct.sasa_1res, ct.rsa_1res, ct.sasa_win3, ct.rsa_win3,
                ct.exposure_class, ct.orientation_angle_deg, ct.total_risk_score, ct.risk_level,
                imgs[0], imgs[1], imgs[2],
            ])


# =============================================================================
# 8. CLI エントリポイント (Main CLI)
# =============================================================================

def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "de novo 酵素と基質・補酵素がドッキングした複数のPDBファイルを一括解析し、"
            "N末・C末のHisタグ付加リスク評価（距離・溶媒露出度RSA）とPyMOL 3方向作図を行いExcelに出力します。"
        )
    )
    parser.add_argument(
        "-i", "--pdb-dir",
        type=Path,
        default=Path("./pdbs"),
        help="解析対象の .pdb ファイルが格納されたディレクトリパス (デフォルト: ./pdbs)",
    )
    parser.add_argument(
        "-o", "--output-excel",
        type=Path,
        default=Path("./histag_risk_evaluation.xlsx"),
        help="出力する Excel (.xlsx) ファイルパス (デフォルト: ./histag_risk_evaluation.xlsx)",
    )
    parser.add_argument(
        "--image-dir",
        type=Path,
        default=Path("./output_images"),
        help="PyMOL で描画した3方向のPNG画像および .pml スクリプトの保存先ディレクトリ (デフォルト: ./output_images)",
    )
    parser.add_argument(
        "--pymol-cmd",
        type=str,
        default=None,
        help="PyMOL 実行ファイルのパス (未指定時は Windows / WSL2 環境から自動検出)",
    )
    parser.add_argument(
        "--substrate-resn",
        type=str,
        nargs="*",
        default=None,
        help="基質 (Substrate) の残基名を明示指定する場合に指定 (例: --substrate-resn SUB LIG UNL)",
    )
    parser.add_argument(
        "--cofactor-resn",
        type=str,
        nargs="*",
        default=None,
        help="補酵素 (Cofactor) の残基名を明示指定する場合に指定 (例: --cofactor-resn NAD FAD HEM SAM)",
    )
    parser.add_argument(
        "--rot-mode",
        choices=["y90_x90", "y90_y180"],
        default="y90_x90",
        help=(
            "追加2種の90度回転図の回転方式: "
            "'y90_x90' = 正面(0°), Y軸90°回転, X軸90°回転の直交3方向 (デフォルト) / "
            "'y90_y180' = 正面(0°), Y軸90°回転, Y軸180°回転"
        ),
    )
    parser.add_argument(
        "--width",
        type=int,
        default=1200,
        help="PyMOL 出力画像の幅 (px, デフォルト: 1200)",
    )
    parser.add_argument(
        "--height",
        type=int,
        default=900,
        help="PyMOL 出力画像の高さ (px, デフォルト: 900)",
    )
    parser.add_argument(
        "--no-ray",
        action="store_true",
        help="高速描画のために PyMOL の ray トレーシングをオフにする (OpenGL 描画)",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    pdb_dir: Path = args.pdb_dir

    if not pdb_dir.exists() or not pdb_dir.is_dir():
        print(f"[ERROR] 指定された PDB ディレクトリが存在しません: {pdb_dir}", file=sys.stderr)
        return 1

    pdb_files = sorted(pdb_dir.glob("*.pdb"))
    if not pdb_files:
        print(f"[ERROR] ディレクトリ内に .pdb ファイルが見つかりません: {pdb_dir}", file=sys.stderr)
        return 1

    pymol_exe = find_pymol_executable(args.pymol_cmd)
    if pymol_exe:
        print(f"[INFO] PyMOL 実行ファイルを検出しました: {pymol_exe}")
    else:
        print(
            "[INFO] PyMOL 実行ファイルが未検出のため、各PDB用の PyMOL スクリプト (.pml) を出力しつつ "
            "内蔵3Dレンダラーで3方向の図を生成します。（PyMOL利用時は --pymol-cmd または環境変数 PYMOL_PATH でも指定可能）"
        )

    results: List[PDBEvaluationResult] = []
    for idx, pdb_path in enumerate(pdb_files, start=1):
        print(f"[{idx}/{len(pdb_files)}] 解析中: {pdb_path.name} ...")
        res = evaluate_pdb_structure(
            pdb_path,
            substrate_resnames=args.substrate_resn,
            cofactor_resnames=args.cofactor_resn,
        )
        render_three_views_for_pdb(
            result=res,
            output_img_dir=args.image_dir,
            pymol_exe=pymol_exe,
            rotation_mode=args.rot_mode,
            image_width=args.width,
            image_height=args.height,
            ray_trace=not args.no_ray,
            allow_fallback=True,
        )
        print(
            f"   -> N末 ({res.n_term.residue_label}): d_min={res.n_term.dist_min_all_ligands} Å, "
            f"RSA_3win={res.n_term.rsa_win3*100:.1f}%, Score={res.n_term.total_risk_score} [{res.n_term.risk_level}]"
        )
        print(
            f"   -> C末 ({res.c_term.residue_label}): d_min={res.c_term.dist_min_all_ligands} Å, "
            f"RSA_3win={res.c_term.rsa_win3*100:.1f}%, Score={res.c_term.total_risk_score} [{res.c_term.risk_level}]"
        )
        print(f"   -> 推奨: {res.recommended_terminus}")
        results.append(res)

    export_to_excel(results, args.output_excel)
    csv_path = args.output_excel.with_suffix(".csv")
    export_to_csv(results, csv_path)

    print(f"\n[DONE] Excel レポートを出力しました: {args.output_excel.resolve()}")
    print(f"[DONE] CSV サマリーを出力しました  : {csv_path.resolve()}")
    print(f"[DONE] 3方向構造画像 & .pml 保存先 : {args.image_dir.resolve()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
