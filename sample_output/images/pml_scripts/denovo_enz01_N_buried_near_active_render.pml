# Auto-generated PyMOL script for denovo_enz01_N_buried_near_active
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

load "C:/Users/a1833/.gemini/antigravity/scratch/histag-risk-evaluator/sample_pdbs/denovo_enz01_N_buried_near_active.pdb", complex
hide everything, all

# 1. セレクション定義
select enzyme, polymer.protein
select ligands, ((chain A and resn NAD and resi 501) or (chain A and resn SUB and resi 502)) and not polymer.protein
select n_term_res, enzyme and ((chain A and resi 1))
select c_term_res, enzyme and ((chain A and resi 40))
select n_term_region, enzyme and ((chain A and resi 1-3))
select c_term_region, enzyme and ((chain A and resi 38-40))

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
png "C:/Users/a1833/.gemini/antigravity/scratch/histag-risk-evaluator/sample_output/images/denovo_enz01_N_buried_near_active_view1_0deg.png", width=1200, height=900, dpi=300, ray=1

# View 2: Rotated 90 deg (#1)
turn y, 90
png "C:/Users/a1833/.gemini/antigravity/scratch/histag-risk-evaluator/sample_output/images/denovo_enz01_N_buried_near_active_view2_rot90_1.png", width=1200, height=900, dpi=300, ray=1

# View 3: Rotated 90 deg (#2)
turn y, -90
turn x, 90
png "C:/Users/a1833/.gemini/antigravity/scratch/histag-risk-evaluator/sample_output/images/denovo_enz01_N_buried_near_active_view3_rot90_2.png", width=1200, height=900, dpi=300, ray=1

quit
